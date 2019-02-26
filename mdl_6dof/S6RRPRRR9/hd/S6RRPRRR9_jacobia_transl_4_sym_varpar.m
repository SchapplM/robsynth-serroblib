% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRR9
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR9_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR9_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_jacobia_transl_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:58:45
% EndTime: 2019-02-26 21:58:45
% DurationCPUTime: 0.10s
% Computational Cost: add. (102->33), mult. (168->53), div. (0->0), fcn. (207->10), ass. (0->28)
t30 = r_i_i_C(3) + pkin(9) + qJ(3);
t14 = sin(pkin(6));
t17 = sin(qJ(2));
t29 = t14 * t17;
t18 = sin(qJ(1));
t28 = t14 * t18;
t20 = cos(qJ(1));
t27 = t14 * t20;
t26 = t17 * t18;
t25 = t17 * t20;
t19 = cos(qJ(2));
t24 = t18 * t19;
t23 = t19 * t20;
t22 = sin(pkin(12)) * pkin(3) + pkin(8);
t12 = pkin(12) + qJ(4);
t10 = sin(t12);
t11 = cos(t12);
t9 = cos(pkin(12)) * pkin(3) + pkin(2);
t21 = t11 * r_i_i_C(1) - t10 * r_i_i_C(2) + t9;
t15 = cos(pkin(6));
t7 = t10 * t27;
t6 = -t15 * t26 + t23;
t5 = t15 * t24 + t25;
t4 = t15 * t25 + t24;
t3 = -t15 * t23 + t26;
t2 = t10 * t28 + t11 * t6;
t1 = -t10 * t6 + t11 * t28;
t8 = [-t18 * pkin(1) + t7 * r_i_i_C(1) - t21 * t4 - t30 * t3 + (r_i_i_C(2) * t11 + t22) * t27, -t21 * t5 + t30 * t6, t5, r_i_i_C(1) * t1 - t2 * r_i_i_C(2), 0, 0; pkin(1) * t20 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t22 * t28 + t30 * t5 + t6 * t9, -t21 * t3 + t30 * t4, t3 (-t4 * t10 - t11 * t27) * r_i_i_C(1) + (-t11 * t4 + t7) * r_i_i_C(2), 0, 0; 0 (t30 * t17 + t21 * t19) * t14, -t14 * t19 (-t10 * t29 + t11 * t15) * r_i_i_C(1) + (-t10 * t15 - t11 * t29) * r_i_i_C(2), 0, 0;];
Ja_transl  = t8;
