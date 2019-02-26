% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRR5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:56:15
% EndTime: 2019-02-26 19:56:15
% DurationCPUTime: 0.18s
% Computational Cost: add. (250->43), mult. (416->74), div. (0->0), fcn. (521->12), ass. (0->35)
t52 = pkin(10) + r_i_i_C(3);
t30 = sin(qJ(6));
t33 = cos(qJ(6));
t51 = t33 * r_i_i_C(1) - t30 * r_i_i_C(2) + pkin(5);
t26 = sin(pkin(11));
t27 = sin(pkin(6));
t48 = t26 * t27;
t28 = cos(pkin(11));
t47 = t27 * t28;
t31 = sin(qJ(4));
t46 = t27 * t31;
t32 = sin(qJ(2));
t45 = t27 * t32;
t35 = cos(qJ(2));
t44 = t27 * t35;
t29 = cos(pkin(6));
t43 = t29 * t32;
t42 = t29 * t35;
t41 = t30 * r_i_i_C(1) + t33 * r_i_i_C(2) + pkin(2) + pkin(8) + pkin(9);
t19 = t26 * t42 + t28 * t32;
t25 = qJ(4) + qJ(5);
t23 = sin(t25);
t24 = cos(t25);
t8 = t19 * t23 + t24 * t48;
t40 = t52 * t8 + t51 * (t19 * t24 - t23 * t48);
t17 = t26 * t32 - t28 * t42;
t10 = -t17 * t23 + t24 * t47;
t39 = -t52 * t10 + t51 * (t17 * t24 + t23 * t47);
t16 = -t23 * t44 + t29 * t24;
t38 = t52 * t16 + t51 * (-t29 * t23 - t24 * t44);
t37 = t31 * pkin(4) + t51 * t23 - t52 * t24 + qJ(3);
t34 = cos(qJ(4));
t20 = -t26 * t43 + t28 * t35;
t18 = t26 * t35 + t28 * t43;
t1 = [0, -t41 * t19 + t37 * t20, t19 (t19 * t34 - t26 * t46) * pkin(4) + t40, t40 (t20 * t33 - t8 * t30) * r_i_i_C(1) + (-t20 * t30 - t8 * t33) * r_i_i_C(2); 0, -t41 * t17 + t37 * t18, t17 (t17 * t34 + t28 * t46) * pkin(4) + t39, t39 (t10 * t30 + t18 * t33) * r_i_i_C(1) + (t10 * t33 - t18 * t30) * r_i_i_C(2); 1 (t37 * t32 + t41 * t35) * t27, -t44 (-t29 * t31 - t34 * t44) * pkin(4) + t38, t38 (-t16 * t30 + t33 * t45) * r_i_i_C(1) + (-t16 * t33 - t30 * t45) * r_i_i_C(2);];
Ja_transl  = t1;
