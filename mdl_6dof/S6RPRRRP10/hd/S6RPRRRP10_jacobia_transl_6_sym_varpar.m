% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP10
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
% Datum: 2019-02-26 21:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRP10_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP10_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:13:04
% EndTime: 2019-02-26 21:13:04
% DurationCPUTime: 0.14s
% Computational Cost: add. (171->34), mult. (205->43), div. (0->0), fcn. (236->8), ass. (0->33)
t19 = sin(qJ(3));
t22 = cos(qJ(3));
t39 = r_i_i_C(2) + pkin(9) + pkin(8);
t21 = cos(qJ(4));
t14 = t21 * pkin(4) + pkin(3);
t17 = qJ(4) + qJ(5);
t15 = sin(t17);
t16 = cos(t17);
t32 = r_i_i_C(3) + qJ(6);
t41 = pkin(5) + r_i_i_C(1);
t42 = t32 * t15 + t41 * t16 + t14;
t45 = t39 * t19 + t42 * t22;
t43 = t39 * t22;
t18 = sin(qJ(4));
t40 = pkin(4) * t18;
t37 = t18 * t19;
t20 = sin(qJ(1));
t36 = t20 * t15;
t35 = t20 * t16;
t23 = cos(qJ(1));
t34 = t23 * t15;
t33 = t23 * t16;
t31 = t32 * t16 * t22;
t30 = t41 * t15;
t29 = pkin(1) + pkin(7) + t40;
t7 = t19 * t36 - t33;
t8 = t19 * t35 + t34;
t28 = t32 * t8 - t41 * t7;
t10 = t19 * t33 - t36;
t9 = t19 * t34 + t35;
t27 = -t32 * t10 + t41 * t9;
t26 = t19 * t14 + qJ(2) - t43;
t1 = [t41 * t10 - t29 * t20 + t26 * t23 + t32 * t9, t20, t45 * t20 (-t20 * t37 + t21 * t23) * pkin(4) + t28, t28, t7; t26 * t20 + t29 * t23 + t32 * t7 + t41 * t8, -t23, -t45 * t23 (t20 * t21 + t23 * t37) * pkin(4) + t27, t27, -t9; 0, 0, -t19 * t42 + t43 (-t30 - t40) * t22 + t31, -t22 * t30 + t31, t22 * t15;];
Ja_transl  = t1;
