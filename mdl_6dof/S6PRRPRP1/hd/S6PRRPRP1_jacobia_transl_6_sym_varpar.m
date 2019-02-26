% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRP1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRP1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:01:10
% EndTime: 2019-02-26 20:01:11
% DurationCPUTime: 0.19s
% Computational Cost: add. (223->43), mult. (368->74), div. (0->0), fcn. (461->12), ass. (0->36)
t46 = pkin(5) + r_i_i_C(1);
t18 = qJ(3) + pkin(11);
t16 = sin(t18);
t17 = cos(t18);
t28 = cos(qJ(3));
t24 = sin(qJ(5));
t27 = cos(qJ(5));
t33 = -t24 * r_i_i_C(2) + t46 * t27 + pkin(4);
t44 = r_i_i_C(3) + qJ(6) + pkin(9);
t45 = t28 * pkin(3) + t44 * t16 + t33 * t17 + pkin(2);
t19 = sin(pkin(10));
t20 = sin(pkin(6));
t43 = t19 * t20;
t21 = cos(pkin(10));
t42 = t20 * t21;
t26 = sin(qJ(2));
t41 = t20 * t26;
t40 = t20 * t28;
t29 = cos(qJ(2));
t39 = t20 * t29;
t38 = cos(pkin(6));
t37 = t26 * t38;
t36 = t29 * t38;
t31 = t27 * r_i_i_C(2) + t46 * t24 + pkin(8) + qJ(4);
t25 = sin(qJ(3));
t10 = -t19 * t37 + t21 * t29;
t9 = t19 * t36 + t21 * t26;
t8 = t19 * t29 + t21 * t37;
t7 = t19 * t26 - t21 * t36;
t6 = t38 * t16 + t17 * t41;
t5 = t16 * t41 - t38 * t17;
t4 = t10 * t17 + t16 * t43;
t3 = t10 * t16 - t17 * t43;
t2 = -t16 * t42 + t8 * t17;
t1 = t8 * t16 + t17 * t42;
t11 = [0, t31 * t10 - t45 * t9, t44 * t4 + (-t10 * t25 + t19 * t40) * pkin(3) - t33 * t3, t9 (-t9 * t24 - t4 * t27) * r_i_i_C(2) + t46 * (-t4 * t24 + t9 * t27) t3; 0, t31 * t8 - t45 * t7, t44 * t2 + (-t21 * t40 - t25 * t8) * pkin(3) - t33 * t1, t7 (-t2 * t27 - t7 * t24) * r_i_i_C(2) + t46 * (-t2 * t24 + t7 * t27) t1; 1 (t31 * t26 + t29 * t45) * t20, t44 * t6 + (-t25 * t41 + t38 * t28) * pkin(3) - t33 * t5, -t39 (t24 * t39 - t6 * t27) * r_i_i_C(2) + t46 * (-t6 * t24 - t27 * t39) t5;];
Ja_transl  = t11;
