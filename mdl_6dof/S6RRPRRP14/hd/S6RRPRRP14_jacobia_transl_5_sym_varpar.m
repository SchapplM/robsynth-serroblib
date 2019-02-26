% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRP14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRP14_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP14_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:53:34
% EndTime: 2019-02-26 21:53:34
% DurationCPUTime: 0.13s
% Computational Cost: add. (164->47), mult. (410->78), div. (0->0), fcn. (520->10), ass. (0->34)
t41 = pkin(2) + pkin(9);
t40 = r_i_i_C(3) + pkin(10);
t19 = sin(pkin(6));
t22 = sin(qJ(2));
t39 = t19 * t22;
t23 = sin(qJ(1));
t38 = t19 * t23;
t26 = cos(qJ(2));
t37 = t19 * t26;
t27 = cos(qJ(1));
t36 = t19 * t27;
t35 = cos(pkin(6));
t34 = t19 * (pkin(3) + pkin(8));
t33 = t23 * t35;
t32 = t27 * t35;
t20 = sin(qJ(5));
t24 = cos(qJ(5));
t31 = t24 * r_i_i_C(1) - t20 * r_i_i_C(2) + pkin(4);
t13 = t23 * t22 - t26 * t32;
t21 = sin(qJ(4));
t25 = cos(qJ(4));
t7 = -t13 * t21 + t25 * t36;
t30 = t13 * t25 + t21 * t36;
t29 = t20 * r_i_i_C(1) + t24 * r_i_i_C(2) + t41;
t28 = t31 * t21 - t40 * t25 + qJ(3);
t16 = -t22 * t33 + t27 * t26;
t15 = t27 * t22 + t26 * t33;
t14 = t22 * t32 + t23 * t26;
t12 = -t21 * t37 + t35 * t25;
t4 = t15 * t21 + t25 * t38;
t3 = -t15 * t25 + t21 * t38;
t2 = t16 * t20 + t4 * t24;
t1 = t16 * t24 - t4 * t20;
t5 = [-t23 * pkin(1) - t13 * qJ(3) - t29 * t14 + t27 * t34 + t40 * t30 + t31 * t7, -t29 * t15 + t28 * t16, t15, -t31 * t3 + t40 * t4, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t27 * pkin(1) + t4 * pkin(4) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t15 * qJ(3) + t41 * t16 + t23 * t34 + t40 * t3, -t29 * t13 + t28 * t14, t13, t31 * t30 - t40 * t7 (t14 * t24 + t7 * t20) * r_i_i_C(1) + (-t14 * t20 + t7 * t24) * r_i_i_C(2), 0; 0 (t28 * t22 + t29 * t26) * t19, -t37, t40 * t12 + t31 * (-t35 * t21 - t25 * t37) (-t12 * t20 + t24 * t39) * r_i_i_C(1) + (-t12 * t24 - t20 * t39) * r_i_i_C(2), 0;];
Ja_transl  = t5;
