% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRP11_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP11_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:14:59
% EndTime: 2019-02-26 22:14:59
% DurationCPUTime: 0.15s
% Computational Cost: add. (183->45), mult. (459->78), div. (0->0), fcn. (583->10), ass. (0->33)
t20 = sin(qJ(3));
t24 = cos(qJ(3));
t19 = sin(qJ(5));
t23 = cos(qJ(5));
t29 = t19 * r_i_i_C(1) + t23 * r_i_i_C(2) + qJ(4);
t33 = pkin(3) + pkin(10) + r_i_i_C(3);
t40 = t29 * t20 + t33 * t24 + pkin(2);
t39 = pkin(4) + pkin(9);
t38 = cos(qJ(1));
t18 = sin(pkin(6));
t22 = sin(qJ(1));
t37 = t18 * t22;
t36 = t18 * t24;
t25 = cos(qJ(2));
t35 = t18 * t25;
t34 = cos(pkin(6));
t32 = t18 * t38;
t31 = t22 * t34;
t30 = t34 * t38;
t28 = t23 * r_i_i_C(1) - t19 * r_i_i_C(2) + t39;
t21 = sin(qJ(2));
t12 = t21 * t30 + t22 * t25;
t3 = t12 * t20 + t24 * t32;
t27 = t12 * t24 - t20 * t32;
t14 = -t21 * t31 + t38 * t25;
t13 = t38 * t21 + t25 * t31;
t11 = t22 * t21 - t25 * t30;
t9 = t18 * t21 * t20 - t34 * t24;
t8 = t14 * t24 + t20 * t37;
t7 = t14 * t20 - t22 * t36;
t2 = t13 * t23 + t7 * t19;
t1 = -t13 * t19 + t7 * t23;
t4 = [-t22 * pkin(1) - t12 * pkin(2) + pkin(8) * t32 - t28 * t11 - t33 * t27 - t29 * t3, -t13 * t40 + t28 * t14, t29 * t8 - t33 * t7, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t38 * pkin(1) + t14 * pkin(2) + pkin(8) * t37 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t7 * qJ(4) + t39 * t13 + t33 * t8, -t11 * t40 + t28 * t12, t29 * t27 - t33 * t3, t3 (-t11 * t19 + t3 * t23) * r_i_i_C(1) + (-t11 * t23 - t3 * t19) * r_i_i_C(2), 0; 0 (t28 * t21 + t40 * t25) * t18, -t33 * t9 + t29 * (t34 * t20 + t21 * t36) t9 (t19 * t35 + t9 * t23) * r_i_i_C(1) + (-t9 * t19 + t23 * t35) * r_i_i_C(2), 0;];
Ja_transl  = t4;
