% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRR2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:35:24
% EndTime: 2019-02-26 20:35:24
% DurationCPUTime: 0.11s
% Computational Cost: add. (192->34), mult. (128->45), div. (0->0), fcn. (142->11), ass. (0->32)
t22 = cos(qJ(5));
t10 = pkin(5) * t22 + pkin(4);
t17 = pkin(11) + qJ(4);
t13 = cos(t17);
t11 = sin(t17);
t35 = r_i_i_C(3) + pkin(9) + pkin(8);
t29 = t35 * t11;
t39 = t29 + t13 * t10 + cos(pkin(11)) * pkin(3) + pkin(2);
t18 = qJ(1) + pkin(10);
t14 = cos(t18);
t19 = qJ(5) + qJ(6);
t16 = cos(t19);
t30 = t14 * t16;
t12 = sin(t18);
t15 = sin(t19);
t34 = t12 * t15;
t5 = t13 * t34 + t30;
t31 = t14 * t15;
t33 = t12 * t16;
t6 = -t13 * t33 + t31;
t38 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t7 = -t13 * t31 + t33;
t8 = t13 * t30 + t34;
t37 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t21 = sin(qJ(5));
t36 = pkin(5) * t21;
t32 = t13 * t21;
t28 = pkin(7) + qJ(3) + t36;
t26 = -r_i_i_C(1) * t15 - r_i_i_C(2) * t16;
t25 = r_i_i_C(1) * t16 - r_i_i_C(2) * t15 + t10;
t24 = -t25 * t11 + t35 * t13;
t1 = [-sin(qJ(1)) * pkin(1) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t28 * t14 - t39 * t12, 0, t12, t24 * t14 (t12 * t22 - t14 * t32) * pkin(5) + t37, t37; cos(qJ(1)) * pkin(1) + t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t28 * t12 + t39 * t14, 0, -t14, t24 * t12 (-t12 * t32 - t14 * t22) * pkin(5) + t38, t38; 0, 1, 0, t25 * t13 + t29 (t26 - t36) * t11, t26 * t11;];
Ja_transl  = t1;
