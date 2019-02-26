% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:28:05
% EndTime: 2019-02-26 21:28:06
% DurationCPUTime: 0.15s
% Computational Cost: add. (229->39), mult. (278->55), div. (0->0), fcn. (337->10), ass. (0->31)
t18 = qJ(2) + pkin(10);
t15 = sin(t18);
t16 = cos(t18);
t38 = pkin(3) + pkin(4);
t43 = cos(qJ(2)) * pkin(2) + t15 * qJ(4) + t38 * t16;
t42 = pkin(1) + t43;
t34 = sin(qJ(5));
t35 = cos(qJ(5));
t8 = t15 * t35 - t16 * t34;
t20 = sin(qJ(6));
t23 = cos(qJ(6));
t27 = t23 * r_i_i_C(1) - t20 * r_i_i_C(2) + pkin(5);
t22 = sin(qJ(1));
t3 = t8 * t22;
t37 = -r_i_i_C(3) - pkin(9);
t7 = t15 * t34 + t16 * t35;
t4 = t7 * t22;
t41 = t27 * t3 - t37 * t4;
t24 = cos(qJ(1));
t29 = t24 * t34;
t30 = t24 * t35;
t5 = -t15 * t29 - t16 * t30;
t6 = -t15 * t30 + t16 * t29;
t40 = -t27 * t6 + t37 * t5;
t39 = -t27 * t7 - t37 * t8;
t36 = -pkin(8) + qJ(3) + pkin(7);
t28 = -t20 * r_i_i_C(1) - t23 * r_i_i_C(2);
t25 = -sin(qJ(2)) * pkin(2) + qJ(4) * t16 - t38 * t15;
t2 = -t22 * t20 - t5 * t23;
t1 = t5 * t20 - t22 * t23;
t9 = [-t37 * t3 - t27 * t4 + (t28 + t36) * t24 - t42 * t22, t25 * t24 - t40, t22, t24 * t15, t40, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); -t5 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t36 * t22 + t42 * t24 - t37 * t6, t25 * t22 - t41, -t24, t22 * t15, t41 (-t4 * t20 + t24 * t23) * r_i_i_C(1) + (-t24 * t20 - t4 * t23) * r_i_i_C(2); 0, -t39 + t43, 0, -t16, t39, t28 * t8;];
Ja_transl  = t9;
