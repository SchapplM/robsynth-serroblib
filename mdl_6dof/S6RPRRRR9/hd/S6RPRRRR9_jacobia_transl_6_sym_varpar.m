% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRR9_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR9_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:19:14
% EndTime: 2019-02-26 21:19:14
% DurationCPUTime: 0.12s
% Computational Cost: add. (196->38), mult. (170->50), div. (0->0), fcn. (188->10), ass. (0->32)
t23 = cos(qJ(3));
t33 = r_i_i_C(3) + pkin(10) + pkin(9) + pkin(8);
t38 = t33 * t23;
t21 = sin(qJ(3));
t20 = qJ(4) + qJ(5);
t17 = qJ(6) + t20;
t13 = sin(t17);
t14 = cos(t17);
t16 = cos(t20);
t11 = pkin(5) * t16 + cos(qJ(4)) * pkin(4);
t9 = pkin(3) + t11;
t26 = r_i_i_C(1) * t14 - r_i_i_C(2) * t13 + t9;
t37 = t33 * t21 + t26 * t23;
t24 = cos(qJ(1));
t22 = sin(qJ(1));
t30 = t22 * t13;
t5 = t14 * t24 - t21 * t30;
t29 = t22 * t14;
t6 = t13 * t24 + t21 * t29;
t36 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t31 = t21 * t24;
t7 = t13 * t31 + t29;
t8 = t14 * t31 - t30;
t35 = t7 * r_i_i_C(1) + t8 * r_i_i_C(2);
t15 = sin(t20);
t34 = pkin(5) * t15;
t32 = t21 * t22;
t10 = t34 + sin(qJ(4)) * pkin(4);
t28 = pkin(1) + pkin(7) + t10;
t27 = -r_i_i_C(1) * t13 - r_i_i_C(2) * t14;
t25 = t21 * t9 + qJ(2) - t38;
t1 = [t8 * r_i_i_C(1) - t7 * r_i_i_C(2) - t28 * t22 + t25 * t24, t22, t37 * t22, -t10 * t32 + t11 * t24 + t36 (-t15 * t32 + t16 * t24) * pkin(5) + t36, t36; t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t25 * t22 + t28 * t24, -t24, -t37 * t24, t10 * t31 + t22 * t11 + t35 (t15 * t31 + t16 * t22) * pkin(5) + t35, t35; 0, 0, -t26 * t21 + t38 (-t10 + t27) * t23 (t27 - t34) * t23, t27 * t23;];
Ja_transl  = t1;
