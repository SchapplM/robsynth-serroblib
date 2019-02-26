% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR6_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR6_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:56:57
% EndTime: 2019-02-26 21:56:57
% DurationCPUTime: 0.13s
% Computational Cost: add. (119->29), mult. (170->41), div. (0->0), fcn. (195->8), ass. (0->28)
t20 = sin(qJ(1));
t17 = qJ(4) + qJ(5);
t16 = cos(t17);
t19 = sin(qJ(2));
t15 = sin(t17);
t22 = cos(qJ(2));
t33 = t15 * t22;
t29 = -t16 * t19 + t33;
t5 = t29 * t20;
t28 = t15 * t19 + t22 * t16;
t6 = t28 * t20;
t37 = -t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t23 = cos(qJ(1));
t32 = t23 * t19;
t7 = -t16 * t32 + t23 * t33;
t8 = t28 * t23;
t36 = -t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t18 = sin(qJ(4));
t30 = pkin(4) * t18 + qJ(3);
t21 = cos(qJ(4));
t35 = pkin(4) * t21 + pkin(2) + pkin(3);
t26 = t19 * t30 + t22 * t35;
t38 = pkin(1) + t26;
t34 = -t28 * r_i_i_C(1) + t29 * r_i_i_C(2);
t31 = pkin(7) - r_i_i_C(3) - pkin(9) - pkin(8);
t27 = pkin(4) * (-t18 * t22 + t19 * t21);
t25 = -t19 * t35 + t22 * t30;
t1 = [-t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t38 * t20 + t31 * t23, t23 * t25 - t36, t32, t23 * t27 + t36, t36, 0; t8 * r_i_i_C(1) - t7 * r_i_i_C(2) + t31 * t20 + t38 * t23, t20 * t25 - t37, t20 * t19, t20 * t27 + t37, t37, 0; 0, t26 - t34, -t22 (-t18 * t19 - t21 * t22) * pkin(4) + t34, t34, 0;];
Ja_transl  = t1;
