% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRR13_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR13_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:55:43
% EndTime: 2019-02-26 20:55:43
% DurationCPUTime: 0.25s
% Computational Cost: add. (453->66), mult. (1247->111), div. (0->0), fcn. (1643->14), ass. (0->52)
t67 = pkin(4) + pkin(9);
t37 = cos(pkin(6));
t35 = cos(pkin(12));
t44 = cos(qJ(1));
t56 = t44 * t35;
t32 = sin(pkin(12));
t41 = sin(qJ(1));
t61 = t41 * t32;
t26 = -t37 * t56 + t61;
t58 = t44 * t32;
t59 = t41 * t35;
t27 = t37 * t58 + t59;
t33 = sin(pkin(7));
t36 = cos(pkin(7));
t40 = sin(qJ(3));
t34 = sin(pkin(6));
t57 = t44 * t34;
t63 = cos(qJ(3));
t66 = (t26 * t36 + t33 * t57) * t40 - t27 * t63;
t39 = sin(qJ(5));
t43 = cos(qJ(5));
t38 = sin(qJ(6));
t42 = cos(qJ(6));
t50 = t42 * r_i_i_C(1) - t38 * r_i_i_C(2) + pkin(5);
t64 = r_i_i_C(3) + pkin(11);
t45 = t50 * t39 - t64 * t43 + qJ(4);
t65 = pkin(10) + pkin(3);
t62 = t34 * t35;
t60 = t41 * t34;
t55 = t34 * qJ(2);
t53 = t33 * t63;
t52 = t36 * t63;
t51 = t34 * t53;
t20 = -t26 * t33 + t36 * t57;
t48 = t37 * t59 + t58;
t47 = t38 * r_i_i_C(1) + t42 * r_i_i_C(2) + t65;
t46 = t48 * t36;
t22 = t33 * t48 + t36 * t60;
t12 = t26 * t52 + t27 * t40 + t44 * t51;
t28 = -t37 * t61 + t56;
t25 = -t33 * t62 + t37 * t36;
t19 = t37 * t33 * t40 + (t35 * t36 * t40 + t63 * t32) * t34;
t18 = t34 * t32 * t40 - t37 * t53 - t52 * t62;
t17 = t28 * t63 + (t33 * t60 - t46) * t40;
t16 = t28 * t40 - t41 * t51 + t63 * t46;
t11 = t18 * t39 + t25 * t43;
t8 = t16 * t39 + t22 * t43;
t7 = -t16 * t43 + t22 * t39;
t6 = t12 * t39 - t20 * t43;
t2 = t17 * t38 + t8 * t42;
t1 = t17 * t42 - t8 * t38;
t3 = [t44 * t55 - t41 * pkin(1) - t27 * pkin(2) + t47 * t66 - t45 * t12 + (t64 * t39 + t50 * t43 + t67) * t20, t60, -t16 * t47 + t17 * t45, t16, -t50 * t7 + t64 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t44 * pkin(1) + t28 * pkin(2) + t8 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t16 * qJ(4) + t65 * t17 + t67 * t22 + t41 * t55 + t64 * t7, -t57, -t12 * t47 - t45 * t66, t12, t64 * t6 + t50 * (t12 * t43 + t20 * t39) (-t6 * t38 - t42 * t66) * r_i_i_C(1) + (t38 * t66 - t6 * t42) * r_i_i_C(2); 0, t37, -t18 * t47 + t19 * t45, t18, t64 * t11 + t50 * (t18 * t43 - t25 * t39) (-t11 * t38 + t19 * t42) * r_i_i_C(1) + (-t11 * t42 - t19 * t38) * r_i_i_C(2);];
Ja_transl  = t3;
