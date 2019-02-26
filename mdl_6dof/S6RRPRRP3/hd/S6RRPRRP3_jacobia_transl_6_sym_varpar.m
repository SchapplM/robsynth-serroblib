% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRP3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:47:09
% EndTime: 2019-02-26 21:47:10
% DurationCPUTime: 0.13s
% Computational Cost: add. (187->39), mult. (157->49), div. (0->0), fcn. (174->10), ass. (0->30)
t21 = qJ(2) + pkin(10);
t14 = sin(t21);
t37 = r_i_i_C(3) + qJ(6) + pkin(9) + pkin(8);
t42 = cos(qJ(2)) * pkin(2) + t37 * t14;
t15 = cos(t21);
t22 = qJ(4) + qJ(5);
t17 = cos(t22);
t11 = pkin(5) * t17 + cos(qJ(4)) * pkin(4);
t9 = pkin(3) + t11;
t41 = t15 * t9 + pkin(1) + t42;
t26 = cos(qJ(1));
t32 = t26 * t17;
t16 = sin(t22);
t25 = sin(qJ(1));
t35 = t25 * t16;
t5 = t15 * t35 + t32;
t33 = t26 * t16;
t34 = t25 * t17;
t6 = -t15 * t34 + t33;
t40 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t7 = -t15 * t33 + t34;
t8 = t15 * t32 + t35;
t39 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t38 = r_i_i_C(2) * t17;
t10 = pkin(5) * t16 + sin(qJ(4)) * pkin(4);
t36 = t10 * t15;
t31 = t10 + qJ(3) + pkin(7);
t28 = r_i_i_C(1) * t17 - r_i_i_C(2) * t16 + t9;
t27 = -sin(qJ(2)) * pkin(2) - t28 * t14 + t37 * t15;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t41 * t25 + t31 * t26, t27 * t26, t25, t25 * t11 - t26 * t36 + t39, t7 * pkin(5) + t39, t26 * t14; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t31 * t25 + t41 * t26, t27 * t25, -t26, -t26 * t11 - t25 * t36 + t40, -t5 * pkin(5) + t40, t25 * t14; 0, t28 * t15 + t42, 0 (-r_i_i_C(1) * t16 - t10 - t38) * t14 (-t38 + (-pkin(5) - r_i_i_C(1)) * t16) * t14, -t15;];
Ja_transl  = t1;
