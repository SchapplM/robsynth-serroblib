% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRR9_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR9_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobia_transl_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:53:23
% EndTime: 2019-02-26 20:53:23
% DurationCPUTime: 0.18s
% Computational Cost: add. (120->56), mult. (315->91), div. (0->0), fcn. (404->12), ass. (0->38)
t26 = sin(qJ(3));
t41 = t26 * pkin(3);
t28 = cos(qJ(3));
t40 = t28 * pkin(3);
t19 = sin(pkin(12));
t27 = sin(qJ(1));
t39 = t27 * t19;
t21 = sin(pkin(6));
t38 = t27 * t21;
t23 = cos(pkin(12));
t37 = t27 * t23;
t29 = cos(qJ(1));
t36 = t29 * t19;
t35 = t29 * t21;
t34 = t29 * t23;
t33 = pkin(9) + qJ(4);
t20 = sin(pkin(7));
t24 = cos(pkin(7));
t18 = sin(pkin(13));
t22 = cos(pkin(13));
t30 = t28 * t18 + t26 * t22;
t4 = t30 * t20;
t32 = r_i_i_C(1) * t4 + t20 * t41 + t33 * t24 + qJ(2);
t25 = cos(pkin(6));
t11 = -t25 * t37 - t36;
t12 = -t25 * t39 + t34;
t13 = t26 * t18 - t28 * t22;
t6 = t30 * t24;
t31 = t11 * t6 - t12 * t13;
t17 = pkin(2) + t40;
t10 = t25 * t36 + t37;
t9 = -t25 * t34 + t39;
t8 = -t33 * t20 + t24 * t41;
t5 = t13 * t24;
t3 = t13 * t20;
t2 = -t11 * t20 + t24 * t38;
t1 = -t11 * t5 - t12 * t30 - t3 * t38;
t7 = [-t27 * pkin(1) + (t13 * r_i_i_C(1) + r_i_i_C(2) * t30 - t17) * t10 + (t6 * r_i_i_C(1) - t5 * r_i_i_C(2) - t20 * r_i_i_C(3) + t8) * t9 + (-r_i_i_C(2) * t3 + r_i_i_C(3) * t24 + t32) * t35, t38, t1 * r_i_i_C(1) + (-t4 * t38 - t31) * r_i_i_C(2) + (-t12 * t26 + (t11 * t24 + t20 * t38) * t28) * pkin(3), t2, 0, 0; t29 * pkin(1) + t31 * r_i_i_C(1) + t1 * r_i_i_C(2) + t2 * r_i_i_C(3) + t11 * t8 + t12 * t17 + t32 * t38, -t35 (-t10 * t30 + t3 * t35 + t9 * t5) * r_i_i_C(1) + (t10 * t13 + t4 * t35 + t9 * t6) * r_i_i_C(2) + (-t10 * t26 + (-t20 * t35 - t24 * t9) * t28) * pkin(3), t9 * t20 - t24 * t35, 0, 0; 0, t25 (-t3 * r_i_i_C(1) - t4 * r_i_i_C(2) + t20 * t40) * t25 + ((-t19 * t30 - t23 * t5) * r_i_i_C(1) + (t13 * t19 - t23 * t6) * r_i_i_C(2) + (t23 * t24 * t28 - t19 * t26) * pkin(3)) * t21, -t21 * t23 * t20 + t25 * t24, 0, 0;];
Ja_transl  = t7;
