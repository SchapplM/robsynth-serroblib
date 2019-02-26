% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function Ja_transl = S6RRRPRP11_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP11_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:14:59
% EndTime: 2019-02-26 22:14:59
% DurationCPUTime: 0.19s
% Computational Cost: add. (233->48), mult. (566->80), div. (0->0), fcn. (716->10), ass. (0->35)
t22 = sin(qJ(3));
t26 = cos(qJ(3));
t21 = sin(qJ(5));
t25 = cos(qJ(5));
t44 = pkin(5) + r_i_i_C(1);
t29 = t25 * r_i_i_C(2) + t21 * t44 + qJ(4);
t37 = pkin(3) + r_i_i_C(3) + qJ(6) + pkin(10);
t28 = -t22 * t29 - t37 * t26 - pkin(2);
t43 = t25 * pkin(5) + pkin(4) + pkin(9);
t42 = cos(qJ(1));
t19 = sin(pkin(6));
t24 = sin(qJ(1));
t41 = t19 * t24;
t40 = t19 * t26;
t27 = cos(qJ(2));
t39 = t19 * t27;
t38 = cos(pkin(6));
t36 = t19 * t42;
t34 = t24 * t38;
t23 = sin(qJ(2));
t13 = t23 * t42 + t27 * t34;
t14 = -t23 * t34 + t27 * t42;
t7 = t14 * t22 - t24 * t40;
t1 = -t13 * t21 + t7 * t25;
t32 = t38 * t42;
t30 = t25 * r_i_i_C(1) - t21 * r_i_i_C(2) + t43;
t12 = t23 * t32 + t24 * t27;
t3 = t12 * t22 + t26 * t36;
t4 = t12 * t26 - t22 * t36;
t11 = t24 * t23 - t27 * t32;
t10 = t22 * t38 + t23 * t40;
t9 = t19 * t23 * t22 - t26 * t38;
t8 = t14 * t26 + t22 * t41;
t2 = t13 * t25 + t7 * t21;
t5 = [-t24 * pkin(1) - t12 * pkin(2) + pkin(8) * t36 - t11 * t30 - t29 * t3 - t37 * t4, t13 * t28 + t14 * t30, t29 * t8 - t37 * t7, t7, -t2 * r_i_i_C(2) + t44 * t1, t8; pkin(8) * t41 + t42 * pkin(1) + t14 * pkin(2) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + (t21 * pkin(5) + qJ(4)) * t7 + t43 * t13 + t37 * t8, t11 * t28 + t12 * t30, t29 * t4 - t3 * t37, t3 (-t11 * t25 - t3 * t21) * r_i_i_C(2) + t44 * (-t11 * t21 + t3 * t25) t4; 0 (t30 * t23 - t28 * t27) * t19, t10 * t29 - t37 * t9, t9 (-t9 * t21 + t25 * t39) * r_i_i_C(2) + t44 * (t21 * t39 + t9 * t25) t10;];
Ja_transl  = t5;
