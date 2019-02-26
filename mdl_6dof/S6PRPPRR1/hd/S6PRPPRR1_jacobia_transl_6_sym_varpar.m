% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPPRR1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPPRR1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:44:35
% EndTime: 2019-02-26 19:44:35
% DurationCPUTime: 0.18s
% Computational Cost: add. (268->46), mult. (540->76), div. (0->0), fcn. (712->13), ass. (0->34)
t26 = sin(pkin(11));
t33 = sin(qJ(2));
t35 = cos(qJ(2));
t43 = cos(pkin(11));
t38 = -t33 * t26 + t35 * t43;
t25 = pkin(12) + qJ(5);
t23 = sin(t25);
t24 = cos(t25);
t32 = sin(qJ(6));
t34 = cos(qJ(6));
t40 = t34 * r_i_i_C(1) - t32 * r_i_i_C(2) + pkin(5);
t48 = pkin(9) + r_i_i_C(3);
t36 = t48 * t23 + t40 * t24 + cos(pkin(12)) * pkin(4) + pkin(3);
t27 = sin(pkin(10));
t28 = sin(pkin(6));
t47 = t27 * t28;
t29 = cos(pkin(10));
t46 = t29 * t28;
t30 = cos(pkin(6));
t45 = t30 * t35;
t19 = -t35 * t26 - t33 * t43;
t17 = t19 * t30;
t7 = -t29 * t17 + t27 * t38;
t41 = -t27 * t17 - t29 * t38;
t39 = t32 * r_i_i_C(1) + t34 * r_i_i_C(2) + pkin(8) + qJ(4);
t37 = t38 * t30;
t16 = t19 * t28;
t15 = t38 * t28;
t12 = -t16 * t24 + t30 * t23;
t9 = t29 * t19 - t27 * t37;
t6 = t27 * t19 + t29 * t37;
t4 = t23 * t47 - t24 * t41;
t2 = -t23 * t46 + t7 * t24;
t1 = [0 (-t27 * t45 - t29 * t33) * pkin(2) - t39 * t41 + t36 * t9, t47, -t9, t48 * t4 + t40 * (t23 * t41 + t24 * t47) (-t4 * t32 - t9 * t34) * r_i_i_C(1) + (t9 * t32 - t4 * t34) * r_i_i_C(2); 0 (-t27 * t33 + t29 * t45) * pkin(2) + t39 * t7 + t36 * t6, -t46, -t6, t48 * t2 + t40 * (-t7 * t23 - t24 * t46) (-t2 * t32 - t6 * t34) * r_i_i_C(1) + (-t2 * t34 + t6 * t32) * r_i_i_C(2); 1, t28 * t35 * pkin(2) + t36 * t15 - t39 * t16, t30, -t15, t48 * t12 + t40 * (t16 * t23 + t30 * t24) (-t12 * t32 - t15 * t34) * r_i_i_C(1) + (-t12 * t34 + t15 * t32) * r_i_i_C(2);];
Ja_transl  = t1;
