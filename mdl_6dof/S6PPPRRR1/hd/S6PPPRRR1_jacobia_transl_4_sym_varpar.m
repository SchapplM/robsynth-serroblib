% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PPPRRR1
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PPPRRR1_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPPRRR1_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_jacobia_transl_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:38:53
% EndTime: 2019-02-26 19:38:53
% DurationCPUTime: 0.13s
% Computational Cost: add. (71->30), mult. (207->65), div. (0->0), fcn. (280->14), ass. (0->37)
t16 = sin(pkin(12));
t19 = sin(pkin(6));
t39 = t16 * t19;
t25 = cos(pkin(6));
t38 = t16 * t25;
t18 = sin(pkin(7));
t37 = t18 * t19;
t36 = t18 * t25;
t21 = cos(pkin(13));
t24 = cos(pkin(7));
t35 = t21 * t24;
t22 = cos(pkin(12));
t34 = t22 * t19;
t33 = t22 * t25;
t15 = sin(pkin(13));
t11 = t15 * t33 + t16 * t21;
t14 = sin(pkin(14));
t17 = sin(pkin(8));
t20 = cos(pkin(14));
t23 = cos(pkin(8));
t10 = -t15 * t16 + t21 * t33;
t29 = t10 * t24 - t18 * t34;
t7 = -t10 * t18 - t24 * t34;
t32 = (-t11 * t14 + t29 * t20) * t23 + t17 * t7;
t13 = -t15 * t38 + t21 * t22;
t12 = -t15 * t22 - t21 * t38;
t28 = t12 * t24 + t16 * t37;
t8 = -t12 * t18 + t24 * t39;
t31 = t17 * t8 + t23 * (-t13 * t14 + t28 * t20);
t9 = -t21 * t37 + t24 * t25;
t30 = t17 * t9 + t23 * (t20 * t36 + (-t14 * t15 + t20 * t35) * t19);
t27 = cos(qJ(4));
t26 = sin(qJ(4));
t6 = t15 * t19 * t20 + (t19 * t35 + t36) * t14;
t4 = t13 * t20 + t28 * t14;
t2 = t11 * t20 + t29 * t14;
t1 = [0, t39, t8 (-t4 * t26 + t31 * t27) * r_i_i_C(1) + (-t31 * t26 - t27 * t4) * r_i_i_C(2), 0, 0; 0, -t34, t7 (-t2 * t26 + t32 * t27) * r_i_i_C(1) + (-t2 * t27 - t32 * t26) * r_i_i_C(2), 0, 0; 1, t25, t9 (-t6 * t26 + t30 * t27) * r_i_i_C(1) + (-t30 * t26 - t27 * t6) * r_i_i_C(2), 0, 0;];
Ja_transl  = t1;
