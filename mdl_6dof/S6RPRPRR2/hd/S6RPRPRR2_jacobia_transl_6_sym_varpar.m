% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRR2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:49:38
% EndTime: 2019-02-26 20:49:38
% DurationCPUTime: 0.12s
% Computational Cost: add. (197->35), mult. (133->46), div. (0->0), fcn. (147->12), ass. (0->32)
t18 = qJ(3) + pkin(11);
t11 = sin(t18);
t37 = r_i_i_C(3) + pkin(9) + pkin(8);
t42 = cos(qJ(3)) * pkin(3) + t37 * t11;
t13 = cos(t18);
t24 = cos(qJ(5));
t9 = pkin(5) * t24 + pkin(4);
t41 = t13 * t9 + pkin(2) + t42;
t19 = qJ(1) + pkin(10);
t14 = cos(t19);
t20 = qJ(5) + qJ(6);
t16 = cos(t20);
t32 = t14 * t16;
t12 = sin(t19);
t15 = sin(t20);
t36 = t12 * t15;
t5 = t13 * t36 + t32;
t33 = t14 * t15;
t35 = t12 * t16;
t6 = -t13 * t35 + t33;
t40 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t7 = -t13 * t33 + t35;
t8 = t13 * t32 + t36;
t39 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t22 = sin(qJ(5));
t38 = pkin(5) * t22;
t34 = t13 * t22;
t30 = qJ(4) + pkin(7) + t38;
t28 = -r_i_i_C(1) * t15 - r_i_i_C(2) * t16;
t27 = r_i_i_C(1) * t16 - r_i_i_C(2) * t15 + t9;
t26 = -sin(qJ(3)) * pkin(3) - t27 * t11 + t37 * t13;
t1 = [-sin(qJ(1)) * pkin(1) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t30 * t14 - t41 * t12, 0, t26 * t14, t12 (t12 * t24 - t14 * t34) * pkin(5) + t39, t39; cos(qJ(1)) * pkin(1) + t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t30 * t12 + t41 * t14, 0, t26 * t12, -t14 (-t12 * t34 - t14 * t24) * pkin(5) + t40, t40; 0, 1, t27 * t13 + t42, 0 (t28 - t38) * t11, t28 * t11;];
Ja_transl  = t1;
