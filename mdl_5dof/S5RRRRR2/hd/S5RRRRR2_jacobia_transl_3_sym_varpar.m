% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S5RRRRR2
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
%
% Output:
% Ja_transl [3x5]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRRRR2_jacobia_transl_3_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_jacobia_transl_3_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR2_jacobia_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_jacobia_transl_3_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-29 15:26:50
% EndTime: 2019-03-29 15:26:51
% DurationCPUTime: 0.06s
% Computational Cost: add. (29->7), mult. (32->12), div. (0->0), fcn. (32->6), ass. (0->10)
t8 = sin(qJ(3));
t9 = cos(qJ(3));
t15 = r_i_i_C(1) * t9 - r_i_i_C(2) * t8;
t12 = -r_i_i_C(1) * t8 - r_i_i_C(2) * t9;
t7 = qJ(1) + qJ(2);
t5 = sin(t7);
t6 = cos(t7);
t11 = t6 * r_i_i_C(3) - t15 * t5;
t10 = t5 * r_i_i_C(3) + t15 * t6;
t1 = [-sin(qJ(1)) * pkin(1) + t11, t11, t12 * t6, 0, 0; cos(qJ(1)) * pkin(1) + t10, t10, t12 * t5, 0, 0; 0, 0, t15, 0, 0;];
Ja_transl  = t1;
