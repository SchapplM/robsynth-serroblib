% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRR1
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
% Datum: 2019-02-26 19:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:53:43
% EndTime: 2019-02-26 19:53:43
% DurationCPUTime: 0.15s
% Computational Cost: add. (148->35), mult. (296->62), div. (0->0), fcn. (381->12), ass. (0->29)
t23 = qJ(4) + qJ(5);
t21 = sin(t23);
t22 = cos(t23);
t26 = sin(pkin(6));
t28 = cos(pkin(11));
t40 = t28 * t26;
t29 = cos(pkin(6));
t24 = sin(pkin(12));
t27 = cos(pkin(12));
t31 = sin(qJ(2));
t33 = cos(qJ(2));
t36 = t33 * t24 + t31 * t27;
t16 = t36 * t29;
t17 = t31 * t24 - t33 * t27;
t25 = sin(pkin(11));
t7 = t28 * t16 - t25 * t17;
t44 = (-t7 * t21 - t22 * t40) * r_i_i_C(1) + (t21 * t40 - t7 * t22) * r_i_i_C(2);
t37 = t25 * t16 + t28 * t17;
t41 = t25 * t26;
t43 = (t21 * t37 + t22 * t41) * r_i_i_C(1) + (-t21 * t41 + t22 * t37) * r_i_i_C(2);
t42 = r_i_i_C(3) + pkin(9) + pkin(8);
t39 = t29 * t33;
t14 = t36 * t26;
t38 = (-t14 * t21 + t29 * t22) * r_i_i_C(1) + (-t14 * t22 - t29 * t21) * r_i_i_C(2);
t32 = cos(qJ(4));
t35 = t32 * pkin(4) + r_i_i_C(1) * t22 - r_i_i_C(2) * t21 + pkin(3);
t30 = sin(qJ(4));
t15 = t17 * t29;
t1 = [0, -t42 * t37 + (-t25 * t39 - t28 * t31) * pkin(2) + t35 * (t25 * t15 - t28 * t36) t41 (t30 * t37 + t32 * t41) * pkin(4) + t43, t43, 0; 0, t42 * t7 + (-t25 * t31 + t28 * t39) * pkin(2) + t35 * (-t28 * t15 - t25 * t36) -t40 (-t30 * t7 - t32 * t40) * pkin(4) + t44, t44, 0; 1, t42 * t14 + (pkin(2) * t33 - t17 * t35) * t26, t29 (-t14 * t30 + t29 * t32) * pkin(4) + t38, t38, 0;];
Ja_transl  = t1;
