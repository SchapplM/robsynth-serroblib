% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRP6_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP6_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:10:49
% EndTime: 2019-02-26 21:10:50
% DurationCPUTime: 0.15s
% Computational Cost: add. (184->38), mult. (152->48), div. (0->0), fcn. (169->9), ass. (0->30)
t20 = pkin(10) + qJ(3);
t15 = cos(t20);
t14 = sin(t20);
t35 = r_i_i_C(3) + qJ(6) + pkin(9) + pkin(8);
t28 = t35 * t14;
t21 = qJ(4) + qJ(5);
t17 = cos(t21);
t11 = pkin(5) * t17 + cos(qJ(4)) * pkin(4);
t9 = pkin(3) + t11;
t39 = t28 + t15 * t9 + cos(pkin(10)) * pkin(2) + pkin(1);
t16 = sin(t21);
t23 = sin(qJ(1));
t31 = t23 * t16;
t24 = cos(qJ(1));
t32 = t17 * t24;
t5 = t15 * t31 + t32;
t30 = t23 * t17;
t33 = t16 * t24;
t6 = -t15 * t30 + t33;
t38 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t7 = -t15 * t33 + t30;
t8 = t15 * t32 + t31;
t37 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t36 = r_i_i_C(2) * t17;
t10 = pkin(5) * t16 + sin(qJ(4)) * pkin(4);
t34 = t10 * t15;
t29 = t10 + pkin(7) + qJ(2);
t26 = r_i_i_C(1) * t17 - r_i_i_C(2) * t16 + t9;
t25 = -t26 * t14 + t35 * t15;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t39 * t23 + t29 * t24, t23, t25 * t24, t23 * t11 - t24 * t34 + t37, t7 * pkin(5) + t37, t24 * t14; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t29 * t23 + t39 * t24, -t24, t25 * t23, -t11 * t24 - t23 * t34 + t38, -t5 * pkin(5) + t38, t23 * t14; 0, 0, t26 * t15 + t28 (-r_i_i_C(1) * t16 - t10 - t36) * t14 (-t36 + (-pkin(5) - r_i_i_C(1)) * t16) * t14, -t15;];
Ja_transl  = t1;
