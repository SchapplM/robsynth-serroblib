% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP9
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRP9_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP9_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:12:30
% EndTime: 2019-02-26 21:12:30
% DurationCPUTime: 0.10s
% Computational Cost: add. (136->37), mult. (152->48), div. (0->0), fcn. (169->8), ass. (0->27)
t20 = cos(qJ(3));
t28 = r_i_i_C(3) + qJ(6) + pkin(9) + pkin(8);
t33 = t28 * t20;
t18 = sin(qJ(3));
t17 = qJ(4) + qJ(5);
t13 = sin(t17);
t14 = cos(t17);
t11 = pkin(5) * t14 + cos(qJ(4)) * pkin(4);
t9 = pkin(3) + t11;
t23 = r_i_i_C(1) * t14 - r_i_i_C(2) * t13 + t9;
t32 = t28 * t18 + t23 * t20;
t21 = cos(qJ(1));
t19 = sin(qJ(1));
t26 = t19 * t13;
t5 = t14 * t21 - t18 * t26;
t25 = t19 * t14;
t6 = t13 * t21 + t18 * t25;
t31 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t27 = t18 * t21;
t7 = t13 * t27 + t25;
t8 = t14 * t27 - t26;
t30 = t7 * r_i_i_C(1) + t8 * r_i_i_C(2);
t29 = r_i_i_C(2) * t14;
t10 = pkin(5) * t13 + sin(qJ(4)) * pkin(4);
t24 = pkin(1) + pkin(7) + t10;
t22 = t18 * t9 + qJ(2) - t33;
t1 = [t8 * r_i_i_C(1) - t7 * r_i_i_C(2) - t24 * t19 + t22 * t21, t19, t32 * t19, -t19 * t18 * t10 + t11 * t21 + t31, t5 * pkin(5) + t31, -t19 * t20; t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t22 * t19 + t24 * t21, -t21, -t32 * t21, t10 * t27 + t19 * t11 + t30, t7 * pkin(5) + t30, t21 * t20; 0, 0, -t23 * t18 + t33 (-r_i_i_C(1) * t13 - t10 - t29) * t20 (-t29 + (-pkin(5) - r_i_i_C(1)) * t13) * t20, t18;];
Ja_transl  = t1;
