% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR4_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR4_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:17:42
% EndTime: 2019-02-26 22:17:42
% DurationCPUTime: 0.12s
% Computational Cost: add. (152->35), mult. (135->47), div. (0->0), fcn. (146->10), ass. (0->29)
t11 = cos(pkin(11)) * pkin(4) + pkin(3);
t20 = qJ(2) + qJ(3);
t16 = sin(t20);
t17 = cos(t20);
t22 = -pkin(9) - qJ(4);
t41 = t17 * t11 + (r_i_i_C(3) - t22) * t16;
t18 = cos(qJ(2)) * pkin(2);
t40 = pkin(1) + t18 + t41;
t24 = sin(qJ(1));
t19 = pkin(11) + qJ(5);
t14 = sin(t19);
t37 = r_i_i_C(2) * t14;
t32 = t16 * t37;
t34 = t17 * t24;
t39 = r_i_i_C(3) * t34 + t24 * t32;
t15 = cos(t19);
t38 = r_i_i_C(1) * t15;
t25 = cos(qJ(1));
t33 = t17 * t25;
t35 = r_i_i_C(3) * t33 + t25 * t32;
t31 = sin(pkin(11)) * pkin(4) + pkin(8) + pkin(7);
t29 = -t17 * t22 + (-t11 - t38) * t16;
t28 = (-t37 + t38) * t17 + t41;
t27 = -sin(qJ(2)) * pkin(2) + t29;
t4 = t14 * t24 + t15 * t33;
t3 = -t14 * t33 + t15 * t24;
t2 = t14 * t25 - t15 * t34;
t1 = t14 * t34 + t15 * t25;
t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t40 * t24 + t31 * t25, t27 * t25 + t35, t29 * t25 + t35, t25 * t16, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t31 * t24 + t40 * t25, t27 * t24 + t39, t29 * t24 + t39, t24 * t16, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0; 0, t18 + t28, t28, -t17 (-r_i_i_C(1) * t14 - r_i_i_C(2) * t15) * t16, 0;];
Ja_transl  = t5;
