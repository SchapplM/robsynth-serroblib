% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRR9_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR9_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobia_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:53:23
% EndTime: 2019-02-26 20:53:23
% DurationCPUTime: 0.22s
% Computational Cost: add. (297->66), mult. (797->111), div. (0->0), fcn. (1047->14), ass. (0->50)
t40 = cos(pkin(6));
t38 = cos(pkin(12));
t46 = cos(qJ(1));
t53 = t46 * t38;
t34 = sin(pkin(12));
t43 = sin(qJ(1));
t58 = t43 * t34;
t24 = -t40 * t53 + t58;
t35 = sin(pkin(7));
t39 = cos(pkin(7));
t36 = sin(pkin(6));
t54 = t46 * t36;
t14 = -t24 * t35 + t39 * t54;
t41 = sin(qJ(5));
t44 = cos(qJ(5));
t33 = sin(pkin(13));
t37 = cos(pkin(13));
t42 = sin(qJ(3));
t45 = cos(qJ(3));
t49 = t45 * t33 + t42 * t37;
t18 = t49 * t35;
t20 = t49 * t39;
t55 = t46 * t34;
t56 = t43 * t38;
t25 = t40 * t55 + t56;
t28 = t42 * t33 - t45 * t37;
t7 = t18 * t54 + t24 * t20 + t25 * t28;
t63 = t14 * t44 - t7 * t41;
t62 = t14 * t41 + t7 * t44;
t13 = t40 * t18 + (t20 * t38 - t28 * t34) * t36;
t61 = -r_i_i_C(3) - pkin(10);
t60 = pkin(3) * t42;
t57 = t43 * t36;
t52 = pkin(9) + qJ(4);
t51 = t36 * (t35 * t60 + t52 * t39 + qJ(2));
t48 = t44 * r_i_i_C(1) - t41 * r_i_i_C(2) + pkin(4);
t17 = t28 * t35;
t19 = t28 * t39;
t47 = -t17 * t54 - t24 * t19 + t25 * t49;
t26 = -t40 * t56 - t55;
t27 = -t40 * t58 + t53;
t10 = t18 * t57 + t26 * t20 - t27 * t28;
t32 = t45 * pkin(3) + pkin(2);
t23 = -t36 * t38 * t35 + t40 * t39;
t22 = -t52 * t35 + t39 * t60;
t16 = -t26 * t35 + t39 * t57;
t9 = -t17 * t57 - t26 * t19 - t27 * t49;
t2 = t10 * t44 + t16 * t41;
t1 = -t10 * t41 + t16 * t44;
t3 = [-t43 * pkin(1) + t7 * pkin(4) + t62 * r_i_i_C(1) + t63 * r_i_i_C(2) + t24 * t22 - t25 * t32 + t46 * t51 + t61 * t47, t57, -t61 * t10 + t48 * t9 + (-t27 * t42 + (t26 * t39 + t35 * t57) * t45) * pkin(3), t16, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t46 * pkin(1) + t10 * pkin(4) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t26 * t22 + t27 * t32 + t43 * t51 + t61 * t9, -t54, t61 * t7 - t48 * t47 + (-t25 * t42 + (-t24 * t39 - t35 * t54) * t45) * pkin(3), -t14, -t63 * r_i_i_C(1) + t62 * r_i_i_C(2), 0; 0, t40, -t61 * t13 + t48 * (-t40 * t17 + (-t19 * t38 - t34 * t49) * t36) + (t35 * t40 * t45 + (t38 * t39 * t45 - t34 * t42) * t36) * pkin(3), t23 (-t13 * t41 + t23 * t44) * r_i_i_C(1) + (-t13 * t44 - t23 * t41) * r_i_i_C(2), 0;];
Ja_transl  = t3;
