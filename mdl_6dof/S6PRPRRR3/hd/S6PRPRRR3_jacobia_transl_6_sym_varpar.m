% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRR3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:54:48
% EndTime: 2019-02-26 19:54:48
% DurationCPUTime: 0.19s
% Computational Cost: add. (345->45), mult. (399->72), div. (0->0), fcn. (500->13), ass. (0->37)
t57 = pkin(10) + r_i_i_C(3);
t33 = sin(qJ(6));
t35 = cos(qJ(6));
t56 = t35 * r_i_i_C(1) - t33 * r_i_i_C(2) + pkin(5);
t29 = pkin(12) + qJ(4);
t27 = qJ(5) + t29;
t23 = sin(t27);
t24 = cos(t27);
t26 = cos(t29);
t55 = t57 * t23 + t56 * t24 + pkin(4) * t26 + cos(pkin(12)) * pkin(3) + pkin(2);
t30 = sin(pkin(11));
t31 = sin(pkin(6));
t51 = t30 * t31;
t34 = sin(qJ(2));
t50 = t30 * t34;
t36 = cos(qJ(2));
t49 = t30 * t36;
t48 = t31 * t34;
t47 = t31 * t36;
t46 = cos(pkin(11));
t45 = t31 * t46;
t44 = t46 * t34;
t43 = t46 * t36;
t41 = t33 * r_i_i_C(1) + t35 * r_i_i_C(2) + pkin(8) + pkin(9) + qJ(3);
t32 = cos(pkin(6));
t17 = t32 * t44 + t49;
t8 = t17 * t24 - t23 * t45;
t40 = t57 * t8 + t56 * (-t17 * t23 - t24 * t45);
t19 = -t32 * t50 + t43;
t10 = t19 * t24 + t23 * t51;
t39 = t57 * t10 + t56 * (-t19 * t23 + t24 * t51);
t15 = t32 * t23 + t24 * t48;
t38 = t57 * t15 + t56 * (-t23 * t48 + t32 * t24);
t25 = sin(t29);
t18 = t32 * t49 + t44;
t16 = -t32 * t43 + t50;
t1 = [0, -t18 * t55 + t41 * t19, t18 (-t19 * t25 + t26 * t51) * pkin(4) + t39, t39 (-t10 * t33 + t18 * t35) * r_i_i_C(1) + (-t10 * t35 - t18 * t33) * r_i_i_C(2); 0, -t16 * t55 + t41 * t17, t16 (-t17 * t25 - t26 * t45) * pkin(4) + t40, t40 (t16 * t35 - t8 * t33) * r_i_i_C(1) + (-t16 * t33 - t8 * t35) * r_i_i_C(2); 1 (t41 * t34 + t55 * t36) * t31, -t47 (-t25 * t48 + t26 * t32) * pkin(4) + t38, t38 (-t15 * t33 - t35 * t47) * r_i_i_C(1) + (-t15 * t35 + t33 * t47) * r_i_i_C(2);];
Ja_transl  = t1;
