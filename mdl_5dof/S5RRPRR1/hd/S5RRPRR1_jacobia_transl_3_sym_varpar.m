% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
%
% Output:
% Ja_transl [3x5]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRPRR1_jacobia_transl_3_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_jacobia_transl_3_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRR1_jacobia_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_jacobia_transl_3_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:22:52
% EndTime: 2019-07-18 17:22:52
% DurationCPUTime: 0.06s
% Computational Cost: add. (15->7), mult. (31->10), div. (0->0), fcn. (33->4), ass. (0->9)
t1 = sin(qJ(2));
t3 = cos(qJ(2));
t8 = pkin(1) + r_i_i_C(1);
t5 = -t1 * r_i_i_C(2) + t8 * t3;
t7 = r_i_i_C(3) + qJ(3);
t6 = -r_i_i_C(2) * t3 - t8 * t1;
t4 = cos(qJ(1));
t2 = sin(qJ(1));
t9 = [-t5 * t2 + t7 * t4, t6 * t4, t2, 0, 0; t7 * t2 + t5 * t4, t6 * t2, -t4, 0, 0; 0, t5, 0, 0, 0;];
Ja_transl  = t9;
