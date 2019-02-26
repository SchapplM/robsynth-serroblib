% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPPRR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPPRR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_jacobia_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:22:34
% EndTime: 2019-02-26 20:22:34
% DurationCPUTime: 0.06s
% Computational Cost: add. (40->10), mult. (30->12), div. (0->0), fcn. (34->6), ass. (0->10)
t9 = -pkin(7) - r_i_i_C(3) + qJ(3);
t4 = sin(qJ(5));
t5 = cos(qJ(5));
t8 = r_i_i_C(1) * t5 - r_i_i_C(2) * t4;
t7 = -t4 * r_i_i_C(1) - t5 * r_i_i_C(2);
t6 = pkin(2) + qJ(4) - t7;
t3 = qJ(1) + pkin(9);
t2 = cos(t3);
t1 = sin(t3);
t10 = [-sin(qJ(1)) * pkin(1) + t9 * t2 - t6 * t1, 0, t1, t2, t8 * t2, 0; cos(qJ(1)) * pkin(1) + t9 * t1 + t6 * t2, 0, -t2, t1, t8 * t1, 0; 0, 1, 0, 0, t7, 0;];
Ja_transl  = t10;
