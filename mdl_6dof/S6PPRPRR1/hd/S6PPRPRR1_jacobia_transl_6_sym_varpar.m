% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PPRPRR1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRPRR1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_jacobia_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:39:47
% EndTime: 2019-02-26 19:39:48
% DurationCPUTime: 0.23s
% Computational Cost: add. (462->59), mult. (1281->112), div. (0->0), fcn. (1706->16), ass. (0->48)
t35 = sin(pkin(7));
t32 = sin(pkin(13));
t37 = cos(pkin(13));
t44 = sin(qJ(3));
t47 = cos(qJ(3));
t52 = t47 * t32 + t44 * t37;
t21 = t52 * t35;
t40 = cos(pkin(7));
t23 = t52 * t40;
t29 = t44 * t32 - t47 * t37;
t33 = sin(pkin(12));
t36 = sin(pkin(6));
t38 = cos(pkin(12));
t41 = cos(pkin(6));
t15 = t41 * t21 + (t23 * t38 - t29 * t33) * t36;
t60 = pkin(10) + r_i_i_C(3);
t34 = sin(pkin(11));
t59 = t34 * t36;
t58 = t34 * t41;
t57 = t35 * t36;
t39 = cos(pkin(11));
t56 = t39 * t36;
t55 = t39 * t41;
t42 = sin(qJ(6));
t45 = cos(qJ(6));
t51 = t45 * r_i_i_C(1) - t42 * r_i_i_C(2) + pkin(5);
t50 = -t42 * r_i_i_C(1) - t45 * r_i_i_C(2) - pkin(9);
t25 = -t34 * t33 + t38 * t55;
t26 = t33 * t55 + t34 * t38;
t49 = t21 * t56 - t25 * t23 + t26 * t29;
t27 = -t39 * t33 - t38 * t58;
t28 = -t33 * t58 + t39 * t38;
t10 = t21 * t59 + t27 * t23 - t28 * t29;
t43 = sin(qJ(5));
t46 = cos(qJ(5));
t48 = t60 * t43 + t51 * t46 + pkin(4);
t24 = -t38 * t57 + t41 * t40;
t22 = t29 * t40;
t20 = t29 * t35;
t17 = -t27 * t35 + t40 * t59;
t16 = -t25 * t35 - t40 * t56;
t14 = -t41 * t20 + (-t22 * t38 - t33 * t52) * t36;
t12 = t15 * t46 + t24 * t43;
t9 = -t20 * t59 - t27 * t22 - t28 * t52;
t6 = t20 * t56 - t25 * t22 - t26 * t52;
t4 = t10 * t46 + t17 * t43;
t2 = t16 * t43 - t46 * t49;
t1 = [0, t59, -t50 * t10 + (-t28 * t44 + (t27 * t40 + t34 * t57) * t47) * pkin(3) + t48 * t9, t17, t60 * t4 + t51 * (-t10 * t43 + t17 * t46) (-t4 * t42 - t9 * t45) * r_i_i_C(1) + (-t4 * t45 + t9 * t42) * r_i_i_C(2); 0, -t56, t50 * t49 + (-t26 * t44 + (t25 * t40 - t35 * t56) * t47) * pkin(3) + t48 * t6, t16, t60 * t2 + t51 * (t16 * t46 + t43 * t49) (-t2 * t42 - t6 * t45) * r_i_i_C(1) + (-t2 * t45 + t6 * t42) * r_i_i_C(2); 1, t41, -t50 * t15 + (t41 * t35 * t47 + (t38 * t40 * t47 - t33 * t44) * t36) * pkin(3) + t48 * t14, t24, t60 * t12 + t51 * (-t15 * t43 + t24 * t46) (-t12 * t42 - t14 * t45) * r_i_i_C(1) + (-t12 * t45 + t14 * t42) * r_i_i_C(2);];
Ja_transl  = t1;
