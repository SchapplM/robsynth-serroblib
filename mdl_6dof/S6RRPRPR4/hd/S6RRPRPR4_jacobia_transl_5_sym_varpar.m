% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR4_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR4_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:39:32
% EndTime: 2019-02-26 21:39:32
% DurationCPUTime: 0.18s
% Computational Cost: add. (186->52), mult. (378->81), div. (0->0), fcn. (487->12), ass. (0->37)
t24 = sin(pkin(11));
t26 = cos(pkin(11));
t30 = sin(qJ(2));
t33 = cos(qJ(2));
t14 = t30 * t24 - t33 * t26;
t29 = sin(qJ(4));
t46 = t29 * pkin(4);
t45 = t33 * pkin(2);
t44 = r_i_i_C(3) + qJ(5) + pkin(9);
t27 = cos(pkin(6));
t43 = t27 * t33;
t25 = sin(pkin(6));
t31 = sin(qJ(1));
t41 = t31 * t25;
t34 = cos(qJ(1));
t39 = t34 * t25;
t37 = t33 * t24 + t30 * t26;
t12 = t37 * t27;
t5 = t34 * t12 - t31 * t14;
t38 = t31 * t12 + t34 * t14;
t32 = cos(qJ(4));
t19 = t32 * pkin(4) + pkin(3);
t23 = qJ(4) + pkin(12);
t21 = sin(t23);
t22 = cos(t23);
t36 = t22 * r_i_i_C(1) - t21 * r_i_i_C(2) + t19;
t35 = t14 * t27;
t20 = pkin(1) + t45;
t16 = t21 * t39;
t13 = t27 * t30 * pkin(2) + (-pkin(8) - qJ(3)) * t25;
t11 = t37 * t25;
t10 = t14 * t25;
t7 = t31 * t35 - t34 * t37;
t4 = -t31 * t37 - t34 * t35;
t2 = t21 * t41 - t22 * t38;
t1 = t21 * t38 + t22 * t41;
t3 = [t16 * r_i_i_C(1) - t31 * t20 - t36 * t5 + t44 * t4 + (-t13 + (t22 * r_i_i_C(2) + t46) * t25) * t34, -t44 * t38 + (-t30 * t34 - t31 * t43) * pkin(2) + t36 * t7, t41, t1 * r_i_i_C(1) - t2 * r_i_i_C(2) + (t29 * t38 + t32 * t41) * pkin(4), -t7, 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t38 * t19 + t34 * t20 - t44 * t7 + (t25 * t46 - t13) * t31, t44 * t5 + (-t30 * t31 + t34 * t43) * pkin(2) + t36 * t4, -t39 (-t5 * t21 - t22 * t39) * r_i_i_C(1) + (-t5 * t22 + t16) * r_i_i_C(2) + (-t5 * t29 - t32 * t39) * pkin(4), -t4, 0; 0, -t10 * t36 + t11 * t44 + t25 * t45, t27 (-t11 * t21 + t27 * t22) * r_i_i_C(1) + (-t11 * t22 - t27 * t21) * r_i_i_C(2) + (-t11 * t29 + t27 * t32) * pkin(4), t10, 0;];
Ja_transl  = t3;
