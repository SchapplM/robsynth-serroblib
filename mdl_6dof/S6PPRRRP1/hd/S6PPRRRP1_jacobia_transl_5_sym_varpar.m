% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PPRRRP1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRRP1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:41:33
% EndTime: 2019-02-26 19:41:33
% DurationCPUTime: 0.17s
% Computational Cost: add. (282->46), mult. (794->87), div. (0->0), fcn. (1045->14), ass. (0->46)
t24 = cos(pkin(11));
t47 = sin(pkin(12));
t48 = sin(pkin(11));
t38 = t48 * t47;
t50 = cos(pkin(12));
t45 = t24 * t50;
t52 = cos(pkin(6));
t32 = -t52 * t45 + t38;
t23 = sin(pkin(6));
t49 = sin(pkin(7));
t46 = t23 * t49;
t51 = cos(pkin(7));
t57 = t24 * t46 + t32 * t51;
t39 = t48 * t50;
t44 = t24 * t47;
t33 = t52 * t39 + t44;
t43 = t48 * t23;
t56 = t33 * t51 - t49 * t43;
t55 = pkin(10) + r_i_i_C(3);
t54 = cos(qJ(3));
t53 = t24 * t23;
t41 = t52 * t49;
t40 = t51 * t50;
t25 = sin(qJ(5));
t28 = cos(qJ(5));
t37 = t28 * r_i_i_C(1) - t25 * r_i_i_C(2) + pkin(4);
t36 = t25 * r_i_i_C(1) + t28 * r_i_i_C(2) + pkin(9);
t26 = sin(qJ(4));
t29 = cos(qJ(4));
t34 = -t55 * t26 - t37 * t29 - pkin(3);
t27 = sin(qJ(3));
t19 = -t52 * t38 + t45;
t18 = t52 * t44 + t39;
t17 = -t50 * t46 + t52 * t51;
t14 = t33 * t49 + t51 * t43;
t13 = t32 * t49 - t51 * t53;
t12 = t27 * t41 + (t27 * t40 + t54 * t47) * t23;
t11 = -t54 * t41 + (t27 * t47 - t40 * t54) * t23;
t10 = t12 * t29 + t17 * t26;
t8 = t19 * t54 - t56 * t27;
t7 = t19 * t27 + t56 * t54;
t6 = t18 * t54 - t57 * t27;
t5 = t18 * t27 + t57 * t54;
t4 = t14 * t26 + t8 * t29;
t2 = t13 * t26 + t6 * t29;
t1 = [0, t43, t34 * t7 + t36 * t8, t55 * t4 + t37 * (t14 * t29 - t8 * t26) (-t4 * t25 + t7 * t28) * r_i_i_C(1) + (-t7 * t25 - t4 * t28) * r_i_i_C(2), 0; 0, -t53, t34 * t5 + t36 * t6, t55 * t2 + t37 * (t13 * t29 - t6 * t26) (-t2 * t25 + t5 * t28) * r_i_i_C(1) + (-t2 * t28 - t5 * t25) * r_i_i_C(2), 0; 1, t52, t34 * t11 + t36 * t12, t55 * t10 + t37 * (-t12 * t26 + t17 * t29) (-t10 * t25 + t11 * t28) * r_i_i_C(1) + (-t10 * t28 - t11 * t25) * r_i_i_C(2), 0;];
Ja_transl  = t1;
