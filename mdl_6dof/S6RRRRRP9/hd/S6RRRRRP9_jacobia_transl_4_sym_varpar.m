% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRRP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRP9_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP9_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_jacobia_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:44:30
% EndTime: 2019-02-26 22:44:30
% DurationCPUTime: 0.17s
% Computational Cost: add. (150->46), mult. (381->80), div. (0->0), fcn. (482->10), ass. (0->34)
t21 = sin(qJ(2));
t22 = sin(qJ(1));
t25 = cos(qJ(2));
t26 = cos(qJ(1));
t33 = cos(pkin(6));
t31 = t26 * t33;
t11 = t22 * t21 - t25 * t31;
t19 = sin(qJ(4));
t23 = cos(qJ(4));
t12 = t21 * t31 + t22 * t25;
t20 = sin(qJ(3));
t24 = cos(qJ(3));
t18 = sin(pkin(6));
t34 = t18 * t26;
t4 = t12 * t24 - t20 * t34;
t43 = -t11 * t23 + t4 * t19;
t42 = -t11 * t19 - t4 * t23;
t30 = t23 * r_i_i_C(1) - t19 * r_i_i_C(2) + pkin(3);
t40 = pkin(10) + r_i_i_C(3);
t41 = t40 * t20 + t30 * t24 + pkin(2);
t37 = t18 * t22;
t36 = t18 * t24;
t35 = t18 * t25;
t32 = t22 * t33;
t29 = t19 * r_i_i_C(1) + t23 * r_i_i_C(2) + pkin(9);
t28 = -t12 * t20 - t24 * t34;
t14 = -t21 * t32 + t26 * t25;
t13 = t26 * t21 + t25 * t32;
t10 = t33 * t20 + t21 * t36;
t8 = t14 * t24 + t20 * t37;
t7 = t14 * t20 - t22 * t36;
t2 = t13 * t19 + t8 * t23;
t1 = t13 * t23 - t8 * t19;
t3 = [-t22 * pkin(1) - t12 * pkin(2) - t4 * pkin(3) + pkin(8) * t34 - t11 * pkin(9) + t42 * r_i_i_C(1) + t43 * r_i_i_C(2) + t40 * t28, -t13 * t41 + t29 * t14, -t30 * t7 + t40 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0; t26 * pkin(1) + t14 * pkin(2) + t8 * pkin(3) + pkin(8) * t37 + t13 * pkin(9) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t40 * t7, -t11 * t41 + t29 * t12, t30 * t28 + t40 * t4, -t43 * r_i_i_C(1) + t42 * r_i_i_C(2), 0, 0; 0 (t29 * t21 + t41 * t25) * t18, t40 * t10 + t30 * (-t18 * t21 * t20 + t33 * t24) (-t10 * t19 - t23 * t35) * r_i_i_C(1) + (-t10 * t23 + t19 * t35) * r_i_i_C(2), 0, 0;];
Ja_transl  = t3;
