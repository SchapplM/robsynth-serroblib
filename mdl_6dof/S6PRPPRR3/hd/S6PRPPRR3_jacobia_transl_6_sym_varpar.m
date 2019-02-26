% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPPRR3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPPRR3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:49
% EndTime: 2019-02-26 19:45:50
% DurationCPUTime: 0.16s
% Computational Cost: add. (212->48), mult. (542->82), div. (0->0), fcn. (715->12), ass. (0->36)
t48 = -pkin(3) - pkin(2);
t47 = pkin(9) + r_i_i_C(3);
t28 = sin(pkin(10));
t29 = sin(pkin(6));
t46 = t28 * t29;
t38 = cos(qJ(2));
t45 = t29 * t38;
t31 = cos(pkin(10));
t44 = t31 * t29;
t32 = cos(pkin(6));
t35 = sin(qJ(2));
t43 = t32 * t35;
t42 = t32 * t38;
t21 = t28 * t35 - t31 * t42;
t22 = t28 * t38 + t31 * t43;
t27 = sin(pkin(11));
t30 = cos(pkin(11));
t7 = -t21 * t30 + t22 * t27;
t8 = t21 * t27 + t22 * t30;
t23 = t28 * t42 + t31 * t35;
t24 = -t28 * t43 + t31 * t38;
t11 = -t23 * t30 + t24 * t27;
t12 = t23 * t27 + t24 * t30;
t33 = sin(qJ(6));
t36 = cos(qJ(6));
t41 = t36 * r_i_i_C(1) - t33 * r_i_i_C(2) + pkin(5);
t40 = t33 * r_i_i_C(1) + t36 * r_i_i_C(2) + pkin(8);
t34 = sin(qJ(5));
t37 = cos(qJ(5));
t39 = t47 * t34 + t41 * t37 + pkin(4);
t20 = t29 * t35 * t30 - t27 * t45;
t19 = (t27 * t35 + t30 * t38) * t29;
t14 = t20 * t37 - t32 * t34;
t4 = t12 * t37 - t34 * t46;
t2 = t34 * t44 + t8 * t37;
t1 = [0, t24 * qJ(3) + t39 * t11 - t40 * t12 + t48 * t23, t23, -t46, t47 * t4 + t41 * (-t12 * t34 - t37 * t46) (t11 * t36 - t4 * t33) * r_i_i_C(1) + (-t11 * t33 - t4 * t36) * r_i_i_C(2); 0, t22 * qJ(3) + t48 * t21 + t39 * t7 - t40 * t8, t21, t44, t47 * t2 + t41 * (-t8 * t34 + t37 * t44) (-t2 * t33 + t7 * t36) * r_i_i_C(1) + (-t2 * t36 - t7 * t33) * r_i_i_C(2); 1 (t35 * qJ(3) - t48 * t38) * t29 - t40 * t20 + t39 * t19, -t45, -t32, t47 * t14 + t41 * (-t20 * t34 - t32 * t37) (-t14 * t33 + t19 * t36) * r_i_i_C(1) + (-t14 * t36 - t19 * t33) * r_i_i_C(2);];
Ja_transl  = t1;
