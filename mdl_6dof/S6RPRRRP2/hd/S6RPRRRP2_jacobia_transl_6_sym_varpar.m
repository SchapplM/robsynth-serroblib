% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP2
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
% Datum: 2019-02-26 21:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRP2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:08:37
% EndTime: 2019-02-26 21:08:37
% DurationCPUTime: 0.10s
% Computational Cost: add. (191->38), mult. (152->51), div. (0->0), fcn. (167->10), ass. (0->28)
t22 = cos(qJ(3));
t21 = sin(qJ(3));
t30 = r_i_i_C(3) + qJ(6) + pkin(9) + pkin(8);
t25 = t30 * t21;
t20 = qJ(4) + qJ(5);
t16 = cos(t20);
t11 = pkin(5) * t16 + cos(qJ(4)) * pkin(4);
t9 = pkin(3) + t11;
t35 = t22 * t9 + pkin(2) + t25;
t19 = qJ(1) + pkin(10);
t13 = sin(t19);
t14 = cos(t19);
t15 = sin(t20);
t28 = t15 * t22;
t5 = t13 * t28 + t14 * t16;
t27 = t16 * t22;
t6 = -t13 * t27 + t14 * t15;
t34 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t7 = t13 * t16 - t14 * t28;
t8 = t13 * t15 + t14 * t27;
t33 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t32 = r_i_i_C(2) * t16;
t10 = pkin(5) * t15 + sin(qJ(4)) * pkin(4);
t31 = pkin(7) + t10;
t29 = t10 * t22;
t24 = r_i_i_C(1) * t16 - r_i_i_C(2) * t15 + t9;
t23 = -t24 * t21 + t30 * t22;
t1 = [-sin(qJ(1)) * pkin(1) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t31 * t14 - t35 * t13, 0, t23 * t14, t13 * t11 - t14 * t29 + t33, t7 * pkin(5) + t33, t14 * t21; cos(qJ(1)) * pkin(1) + t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t31 * t13 + t35 * t14, 0, t23 * t13, -t14 * t11 - t13 * t29 + t34, -t5 * pkin(5) + t34, t13 * t21; 0, 1, t24 * t22 + t25 (-r_i_i_C(1) * t15 - t10 - t32) * t21 (-t32 + (-pkin(5) - r_i_i_C(1)) * t15) * t21, -t22;];
Ja_transl  = t1;
