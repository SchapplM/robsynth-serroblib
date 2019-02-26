% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR2_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR2_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:38:16
% EndTime: 2019-02-26 21:38:16
% DurationCPUTime: 0.09s
% Computational Cost: add. (121->21), mult. (76->22), div. (0->0), fcn. (81->8), ass. (0->18)
t36 = r_i_i_C(3) + qJ(5);
t18 = qJ(2) + pkin(10);
t15 = qJ(4) + t18;
t13 = sin(t15);
t14 = cos(t15);
t21 = (pkin(4) - r_i_i_C(2)) * t14 + t36 * t13;
t35 = t21 + pkin(3) * cos(t18) + cos(qJ(2)) * pkin(2);
t34 = t36 * t14;
t33 = pkin(1) + t35;
t30 = r_i_i_C(1) + pkin(8) + qJ(3) + pkin(7);
t19 = sin(qJ(1));
t29 = t19 * t13;
t20 = cos(qJ(1));
t28 = t20 * t13;
t24 = r_i_i_C(2) * t29 + t34 * t19;
t23 = r_i_i_C(2) * t28 + t34 * t20;
t22 = -pkin(4) * t13 - pkin(3) * sin(t18) - sin(qJ(2)) * pkin(2);
t1 = [-t33 * t19 + t30 * t20, t22 * t20 + t23, t19, -pkin(4) * t28 + t23, t28, 0; t30 * t19 + t33 * t20, t22 * t19 + t24, -t20, -pkin(4) * t29 + t24, t29, 0; 0, t35, 0, t21, -t14, 0;];
Ja_transl  = t1;
