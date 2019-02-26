% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRP5_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP5_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:41:48
% EndTime: 2019-02-26 22:41:48
% DurationCPUTime: 0.12s
% Computational Cost: add. (193->36), mult. (168->50), div. (0->0), fcn. (184->10), ass. (0->32)
t23 = cos(qJ(2));
t21 = sin(qJ(2));
t34 = r_i_i_C(3) + pkin(10) + pkin(9) + pkin(8);
t28 = t34 * t21;
t20 = qJ(3) + qJ(4);
t16 = cos(t20);
t11 = pkin(4) * t16 + cos(qJ(3)) * pkin(3);
t9 = pkin(2) + t11;
t39 = t23 * t9 + pkin(1) + t28;
t17 = qJ(5) + t20;
t13 = sin(t17);
t14 = cos(t17);
t24 = cos(qJ(1));
t30 = t24 * t14;
t22 = sin(qJ(1));
t33 = t22 * t23;
t5 = t13 * t33 + t30;
t31 = t24 * t13;
t6 = -t14 * t33 + t31;
t38 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t7 = t22 * t14 - t23 * t31;
t8 = t22 * t13 + t23 * t30;
t37 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t15 = sin(t20);
t36 = pkin(4) * t15;
t10 = t36 + sin(qJ(3)) * pkin(3);
t35 = pkin(7) + t10;
t32 = t23 * t24;
t27 = -r_i_i_C(1) * t13 - r_i_i_C(2) * t14;
t26 = r_i_i_C(1) * t14 - r_i_i_C(2) * t13 + t9;
t25 = -t26 * t21 + t34 * t23;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t39 * t22 + t35 * t24, t25 * t24, -t10 * t32 + t22 * t11 + t37 (-t15 * t32 + t16 * t22) * pkin(4) + t37, t37, 0; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t35 * t22 + t39 * t24, t25 * t22, -t10 * t33 - t24 * t11 + t38 (-t15 * t33 - t16 * t24) * pkin(4) + t38, t38, 0; 0, t26 * t23 + t28 (-t10 + t27) * t21 (t27 - t36) * t21, t27 * t21, 0;];
Ja_transl  = t1;
