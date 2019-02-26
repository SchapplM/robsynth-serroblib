% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPR8_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR8_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:05:01
% EndTime: 2019-02-26 21:05:01
% DurationCPUTime: 0.10s
% Computational Cost: add. (160->34), mult. (140->45), div. (0->0), fcn. (157->10), ass. (0->28)
t21 = cos(qJ(3));
t30 = r_i_i_C(3) + pkin(9) + qJ(5) + pkin(8);
t34 = t30 * t21;
t19 = sin(qJ(3));
t18 = qJ(4) + pkin(10);
t15 = qJ(6) + t18;
t13 = sin(t15);
t14 = cos(t15);
t11 = pkin(5) * cos(t18) + cos(qJ(4)) * pkin(4);
t9 = pkin(3) + t11;
t24 = r_i_i_C(1) * t14 - r_i_i_C(2) * t13 + t9;
t33 = t30 * t19 + t24 * t21;
t22 = cos(qJ(1));
t20 = sin(qJ(1));
t28 = t20 * t13;
t5 = t14 * t22 - t19 * t28;
t27 = t20 * t14;
t6 = t13 * t22 + t19 * t27;
t32 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t29 = t19 * t22;
t7 = t13 * t29 + t27;
t8 = t14 * t29 - t28;
t31 = t7 * r_i_i_C(1) + t8 * r_i_i_C(2);
t10 = pkin(5) * sin(t18) + sin(qJ(4)) * pkin(4);
t26 = pkin(1) + pkin(7) + t10;
t25 = -r_i_i_C(1) * t13 - r_i_i_C(2) * t14;
t23 = t19 * t9 + qJ(2) - t34;
t1 = [t8 * r_i_i_C(1) - t7 * r_i_i_C(2) - t26 * t20 + t23 * t22, t20, t33 * t20, -t20 * t19 * t10 + t11 * t22 + t32, -t20 * t21, t32; t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t23 * t20 + t26 * t22, -t22, -t33 * t22, t10 * t29 + t20 * t11 + t31, t22 * t21, t31; 0, 0, -t24 * t19 + t34 (-t10 + t25) * t21, t19, t25 * t21;];
Ja_transl  = t1;
