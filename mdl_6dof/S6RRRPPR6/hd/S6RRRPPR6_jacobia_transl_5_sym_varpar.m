% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPR6_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR6_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:06:32
% EndTime: 2019-02-26 22:06:33
% DurationCPUTime: 0.13s
% Computational Cost: add. (185->41), mult. (294->63), div. (0->0), fcn. (365->10), ass. (0->34)
t25 = cos(qJ(3));
t15 = t25 * pkin(3) + pkin(2);
t18 = qJ(3) + pkin(11);
t16 = sin(t18);
t17 = cos(t18);
t31 = r_i_i_C(3) + qJ(5);
t40 = pkin(4) - r_i_i_C(2);
t41 = t31 * t16 + t40 * t17 + t15;
t39 = r_i_i_C(1) + qJ(4) + pkin(9);
t19 = sin(pkin(6));
t23 = sin(qJ(2));
t38 = t19 * t23;
t24 = sin(qJ(1));
t37 = t19 * t24;
t27 = cos(qJ(1));
t36 = t19 * t27;
t35 = t24 * t23;
t26 = cos(qJ(2));
t34 = t24 * t26;
t33 = t27 * t23;
t32 = t27 * t26;
t22 = sin(qJ(3));
t30 = t19 * (pkin(3) * t22 + pkin(8));
t20 = cos(pkin(6));
t10 = t20 * t33 + t34;
t1 = t10 * t16 + t17 * t36;
t29 = -t10 * t17 + t16 * t36;
t12 = -t20 * t35 + t32;
t11 = t20 * t34 + t33;
t9 = -t20 * t32 + t35;
t7 = t16 * t38 - t20 * t17;
t6 = t12 * t17 + t16 * t37;
t5 = t12 * t16 - t17 * t37;
t2 = [-t24 * pkin(1) - t31 * t1 - t10 * t15 + t27 * t30 + t40 * t29 - t39 * t9, -t11 * t41 + t39 * t12, t31 * t6 - t40 * t5 + (-t12 * t22 + t25 * t37) * pkin(3), t11, t5, 0; t27 * pkin(1) + t39 * t11 + t12 * t15 + t24 * t30 + t31 * t5 + t40 * t6, t39 * t10 - t41 * t9, -t31 * t29 - t40 * t1 + (-t10 * t22 - t25 * t36) * pkin(3), t9, t1, 0; 0 (t39 * t23 + t41 * t26) * t19, t31 * (t20 * t16 + t17 * t38) - t40 * t7 + (t20 * t25 - t22 * t38) * pkin(3), -t19 * t26, t7, 0;];
Ja_transl  = t2;
