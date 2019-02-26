% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Ja_transl = S6PPPRRR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPPRRR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_jacobia_transl_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:38:43
% EndTime: 2019-02-26 19:38:44
% DurationCPUTime: 0.18s
% Computational Cost: add. (280->45), mult. (804->91), div. (0->0), fcn. (1074->16), ass. (0->50)
t52 = pkin(10) + r_i_i_C(3);
t25 = sin(pkin(12));
t28 = sin(pkin(6));
t51 = t25 * t28;
t34 = cos(pkin(6));
t50 = t25 * t34;
t27 = sin(pkin(7));
t49 = t27 * t28;
t48 = t27 * t34;
t30 = cos(pkin(13));
t33 = cos(pkin(7));
t47 = t30 * t33;
t31 = cos(pkin(12));
t46 = t31 * t28;
t45 = t31 * t34;
t24 = sin(pkin(13));
t19 = -t25 * t24 + t30 * t45;
t16 = -t19 * t27 - t33 * t46;
t26 = sin(pkin(8));
t32 = cos(pkin(8));
t20 = t24 * t45 + t25 * t30;
t23 = sin(pkin(14));
t29 = cos(pkin(14));
t40 = t19 * t33 - t27 * t46;
t9 = -t20 * t23 + t40 * t29;
t44 = t16 * t26 + t32 * t9;
t22 = -t24 * t50 + t31 * t30;
t21 = -t31 * t24 - t30 * t50;
t39 = t21 * t33 + t25 * t49;
t11 = -t22 * t23 + t39 * t29;
t17 = -t21 * t27 + t33 * t51;
t43 = t11 * t32 + t17 * t26;
t14 = t29 * t48 + (-t23 * t24 + t29 * t47) * t28;
t18 = -t30 * t49 + t34 * t33;
t42 = t14 * t32 + t18 * t26;
t35 = sin(qJ(5));
t37 = cos(qJ(5));
t41 = t37 * r_i_i_C(1) - t35 * r_i_i_C(2) + pkin(4);
t38 = cos(qJ(4));
t36 = sin(qJ(4));
t15 = t28 * t24 * t29 + (t28 * t47 + t48) * t23;
t13 = -t14 * t26 + t18 * t32;
t12 = t22 * t29 + t39 * t23;
t10 = t20 * t29 + t40 * t23;
t8 = -t11 * t26 + t17 * t32;
t7 = t16 * t32 - t9 * t26;
t6 = t15 * t38 + t42 * t36;
t4 = t12 * t38 + t43 * t36;
t2 = t10 * t38 + t44 * t36;
t1 = [0, t51, t17, t52 * t4 + t41 * (-t12 * t36 + t43 * t38) (-t4 * t35 + t8 * t37) * r_i_i_C(1) + (-t8 * t35 - t4 * t37) * r_i_i_C(2), 0; 0, -t46, t16, t52 * t2 + t41 * (-t10 * t36 + t44 * t38) (-t2 * t35 + t7 * t37) * r_i_i_C(1) + (-t2 * t37 - t7 * t35) * r_i_i_C(2), 0; 1, t34, t18, t52 * t6 + t41 * (-t15 * t36 + t42 * t38) (t13 * t37 - t6 * t35) * r_i_i_C(1) + (-t13 * t35 - t6 * t37) * r_i_i_C(2), 0;];
Ja_transl  = t1;
