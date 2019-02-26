% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR3
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRR3_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR3_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_jacobia_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:23
% EndTime: 2019-02-26 20:05:24
% DurationCPUTime: 0.25s
% Computational Cost: add. (288->73), mult. (781->137), div. (0->0), fcn. (1021->14), ass. (0->45)
t33 = sin(pkin(13));
t41 = sin(qJ(3));
t44 = cos(qJ(3));
t50 = cos(pkin(13));
t29 = -t44 * t33 - t41 * t50;
t35 = sin(pkin(7));
t19 = t29 * t35;
t38 = cos(pkin(7));
t21 = t29 * t38;
t36 = sin(pkin(6));
t39 = cos(pkin(6));
t42 = sin(qJ(2));
t45 = cos(qJ(2));
t47 = -t41 * t33 + t44 * t50;
t9 = -t39 * t19 + (-t21 * t45 + t42 * t47) * t36;
t60 = -r_i_i_C(3) - pkin(10);
t34 = sin(pkin(12));
t59 = t34 * t36;
t58 = t35 * t36;
t40 = sin(qJ(5));
t57 = t35 * t40;
t43 = cos(qJ(5));
t56 = t35 * t43;
t37 = cos(pkin(12));
t55 = t36 * t37;
t54 = t36 * t38;
t52 = t39 * t42;
t51 = t39 * t45;
t48 = t43 * r_i_i_C(1) - t40 * r_i_i_C(2) + pkin(4);
t24 = -t34 * t42 + t37 * t51;
t25 = t34 * t45 + t37 * t52;
t46 = -t19 * t55 + t24 * t21 - t25 * t47;
t26 = -t34 * t51 - t37 * t42;
t27 = -t34 * t52 + t37 * t45;
t6 = -t19 * t59 - t26 * t21 + t27 * t47;
t32 = t44 * pkin(3) + pkin(2);
t23 = t39 * t38 - t45 * t58;
t22 = t38 * t41 * pkin(3) + (-pkin(9) - qJ(4)) * t35;
t20 = t47 * t38;
t18 = t47 * t35;
t17 = -t26 * t35 + t34 * t54;
t16 = -t24 * t35 - t37 * t54;
t13 = t27 * t21 + t26 * t47;
t11 = t25 * t21 + t24 * t47;
t1 = [0 (t13 * t43 + t27 * t57) * r_i_i_C(1) + (-t13 * t40 + t27 * t56) * r_i_i_C(2) + t13 * pkin(4) + t26 * t32 - t27 * t22 + t60 * (-t27 * t20 + t26 * t29) -t60 * t6 + t48 * (t18 * t59 + t26 * t20 + t27 * t29) + (-t27 * t41 + (t26 * t38 + t34 * t58) * t44) * pkin(3), t17 (t17 * t43 - t6 * t40) * r_i_i_C(1) + (-t17 * t40 - t6 * t43) * r_i_i_C(2), 0; 0 (t11 * t43 + t25 * t57) * r_i_i_C(1) + (-t11 * t40 + t25 * t56) * r_i_i_C(2) + t11 * pkin(4) + t24 * t32 - t25 * t22 + t60 * (-t25 * t20 + t24 * t29) t60 * t46 + t48 * (-t18 * t55 + t24 * t20 + t25 * t29) + (-t25 * t41 + (t24 * t38 - t35 * t55) * t44) * pkin(3), t16 (t16 * t43 + t40 * t46) * r_i_i_C(1) + (-t16 * t40 + t43 * t46) * r_i_i_C(2), 0; 1 (t48 * (t21 * t42 + t45 * t47) + t45 * t32 + (-t22 + (t40 * r_i_i_C(1) + t43 * r_i_i_C(2)) * t35) * t42 + t60 * (-t20 * t42 + t29 * t45)) * t36, -t60 * t9 + t48 * (t39 * t18 + (t20 * t45 + t29 * t42) * t36) + (t35 * t39 * t44 + (t38 * t44 * t45 - t41 * t42) * t36) * pkin(3), t23 (t23 * t43 - t9 * t40) * r_i_i_C(1) + (-t23 * t40 - t9 * t43) * r_i_i_C(2), 0;];
Ja_transl  = t1;
