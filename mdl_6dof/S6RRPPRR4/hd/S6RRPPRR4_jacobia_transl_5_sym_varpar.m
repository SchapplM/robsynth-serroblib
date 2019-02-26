% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR4_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR4_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:10
% EndTime: 2019-02-26 21:30:11
% DurationCPUTime: 0.16s
% Computational Cost: add. (155->44), mult. (378->69), div. (0->0), fcn. (493->10), ass. (0->31)
t21 = sin(pkin(11));
t23 = cos(pkin(11));
t26 = sin(qJ(2));
t29 = cos(qJ(2));
t15 = t26 * t21 - t29 * t23;
t42 = t29 * pkin(2);
t24 = cos(pkin(6));
t41 = t24 * t29;
t22 = sin(pkin(6));
t27 = sin(qJ(1));
t39 = t27 * t22;
t30 = cos(qJ(1));
t37 = t30 * t22;
t36 = -r_i_i_C(3) - pkin(9) - pkin(3);
t33 = t29 * t21 + t26 * t23;
t13 = t33 * t24;
t35 = t30 * t13 - t27 * t15;
t34 = t27 * t13 + t30 * t15;
t25 = sin(qJ(5));
t28 = cos(qJ(5));
t32 = t25 * r_i_i_C(1) + t28 * r_i_i_C(2) + qJ(4);
t31 = t15 * t24;
t20 = pkin(1) + t42;
t18 = t28 * t37;
t14 = t24 * t26 * pkin(2) + (-pkin(8) - qJ(3)) * t22;
t11 = t15 * t22;
t7 = t27 * t31 - t30 * t33;
t4 = -t27 * t33 - t30 * t31;
t2 = -t7 * t25 + t28 * t39;
t1 = -t25 * t39 - t7 * t28;
t3 = [t18 * r_i_i_C(1) - t27 * t20 + t32 * t4 + (-t14 + (-t25 * r_i_i_C(2) + pkin(4)) * t22) * t30 + t36 * t35 (-t26 * t30 - t27 * t41) * pkin(2) - t36 * t7 - t32 * t34, t39, -t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t7 * qJ(4) + t30 * t20 + (t22 * pkin(4) - t14) * t27 + t36 * t34 (-t26 * t27 + t30 * t41) * pkin(2) - t36 * t4 + t32 * t35, -t37, -t4 (t25 * t37 - t4 * t28) * r_i_i_C(1) + (t4 * t25 + t18) * r_i_i_C(2), 0; 0, t36 * t11 + (t32 * t33 + t42) * t22, t24, t11 (t11 * t28 - t24 * t25) * r_i_i_C(1) + (-t11 * t25 - t24 * t28) * r_i_i_C(2), 0;];
Ja_transl  = t3;
