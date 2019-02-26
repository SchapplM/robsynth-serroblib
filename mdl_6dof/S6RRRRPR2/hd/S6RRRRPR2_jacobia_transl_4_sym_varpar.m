% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR2_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR2_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_jacobia_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:31:15
% EndTime: 2019-02-26 22:31:15
% DurationCPUTime: 0.07s
% Computational Cost: add. (86->12), mult. (61->18), div. (0->0), fcn. (61->8), ass. (0->15)
t11 = qJ(2) + qJ(3);
t8 = qJ(4) + t11;
t5 = sin(t8);
t6 = cos(t8);
t19 = t6 * r_i_i_C(1) - t5 * r_i_i_C(2);
t18 = t19 + pkin(3) * cos(t11);
t23 = t18 + cos(qJ(2)) * pkin(2);
t17 = -r_i_i_C(1) * t5 - r_i_i_C(2) * t6;
t14 = t17 - pkin(3) * sin(t11);
t20 = r_i_i_C(3) + pkin(9) + pkin(8) + pkin(7);
t16 = pkin(1) + t23;
t15 = -sin(qJ(2)) * pkin(2) + t14;
t13 = cos(qJ(1));
t12 = sin(qJ(1));
t1 = [-t16 * t12 + t20 * t13, t15 * t13, t14 * t13, t17 * t13, 0, 0; t20 * t12 + t16 * t13, t15 * t12, t14 * t12, t17 * t12, 0, 0; 0, t23, t18, t19, 0, 0;];
Ja_transl  = t1;
