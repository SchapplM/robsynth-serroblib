% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PPRRRR3_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRRR3_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobia_transl_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:43:55
% EndTime: 2019-02-26 19:43:55
% DurationCPUTime: 0.22s
% Computational Cost: add. (454->73), mult. (1310->140), div. (0->0), fcn. (1733->16), ass. (0->58)
t64 = r_i_i_C(3) + pkin(11);
t32 = sin(pkin(8));
t63 = t32 * pkin(10);
t31 = sin(pkin(13));
t34 = sin(pkin(6));
t62 = t31 * t34;
t39 = cos(pkin(6));
t61 = t31 * t39;
t40 = sin(qJ(5));
t60 = t32 * t40;
t43 = cos(qJ(5));
t59 = t32 * t43;
t33 = sin(pkin(7));
t58 = t33 * t34;
t57 = t33 * t39;
t36 = cos(pkin(13));
t56 = t36 * t34;
t55 = t36 * t39;
t37 = cos(pkin(8));
t41 = sin(qJ(4));
t54 = t37 * t41;
t44 = cos(qJ(4));
t53 = t37 * t44;
t38 = cos(pkin(7));
t42 = sin(qJ(3));
t52 = t38 * t42;
t51 = t33 * t56;
t30 = sin(pkin(14));
t35 = cos(pkin(14));
t25 = -t31 * t30 + t35 * t55;
t26 = t30 * t55 + t31 * t35;
t45 = cos(qJ(3));
t15 = -t26 * t42 + (t25 * t38 - t51) * t45;
t22 = -t25 * t33 - t38 * t56;
t50 = t15 * t37 + t22 * t32;
t28 = -t30 * t61 + t36 * t35;
t27 = -t36 * t30 - t35 * t61;
t46 = t27 * t38 + t31 * t58;
t17 = -t28 * t42 + t46 * t45;
t23 = -t27 * t33 + t38 * t62;
t49 = t17 * t37 + t23 * t32;
t20 = t45 * t57 + (t35 * t38 * t45 - t30 * t42) * t34;
t24 = -t35 * t58 + t39 * t38;
t48 = t20 * t37 + t24 * t32;
t47 = t43 * r_i_i_C(1) - t40 * r_i_i_C(2) + pkin(4);
t21 = t42 * t57 + (t30 * t45 + t35 * t52) * t34;
t19 = -t20 * t32 + t24 * t37;
t18 = t28 * t45 + t46 * t42;
t16 = t25 * t52 + t26 * t45 - t42 * t51;
t14 = t20 * t44 - t21 * t54;
t12 = -t17 * t32 + t23 * t37;
t11 = -t15 * t32 + t22 * t37;
t10 = t21 * t44 + t48 * t41;
t8 = t17 * t44 - t18 * t54;
t6 = t15 * t44 - t16 * t54;
t4 = t18 * t44 + t49 * t41;
t2 = t16 * t44 + t50 * t41;
t1 = [0, t62 (t18 * t60 + t8 * t43) * r_i_i_C(1) + (t18 * t59 - t8 * t40) * r_i_i_C(2) + t8 * pkin(4) + t17 * pkin(3) + t18 * t63 + t64 * (t17 * t41 + t18 * t53) t64 * t4 + t47 * (-t18 * t41 + t49 * t44) (t12 * t43 - t4 * t40) * r_i_i_C(1) + (-t12 * t40 - t4 * t43) * r_i_i_C(2), 0; 0, -t56 (t16 * t60 + t6 * t43) * r_i_i_C(1) + (t16 * t59 - t6 * t40) * r_i_i_C(2) + t6 * pkin(4) + t15 * pkin(3) + t16 * t63 + t64 * (t15 * t41 + t16 * t53) t64 * t2 + t47 * (-t16 * t41 + t50 * t44) (t11 * t43 - t2 * t40) * r_i_i_C(1) + (-t11 * t40 - t2 * t43) * r_i_i_C(2), 0; 1, t39 (t14 * t43 + t21 * t60) * r_i_i_C(1) + (-t14 * t40 + t21 * t59) * r_i_i_C(2) + t14 * pkin(4) + t20 * pkin(3) + t21 * t63 + t64 * (t20 * t41 + t21 * t53) t64 * t10 + t47 * (-t21 * t41 + t48 * t44) (-t10 * t40 + t19 * t43) * r_i_i_C(1) + (-t10 * t43 - t19 * t40) * r_i_i_C(2), 0;];
Ja_transl  = t1;
