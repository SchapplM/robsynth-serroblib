% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR11_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR11_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:21:49
% EndTime: 2019-02-26 22:21:49
% DurationCPUTime: 0.16s
% Computational Cost: add. (199->45), mult. (502->77), div. (0->0), fcn. (645->10), ass. (0->33)
t22 = sin(qJ(3));
t26 = cos(qJ(3));
t21 = sin(qJ(5));
t25 = cos(qJ(5));
t39 = pkin(3) + pkin(4);
t29 = t25 * r_i_i_C(1) - t21 * r_i_i_C(2) + t39;
t30 = t21 * r_i_i_C(1) + t25 * r_i_i_C(2) + qJ(4);
t40 = t30 * t22 + t29 * t26 + pkin(2);
t38 = cos(qJ(1));
t20 = sin(pkin(6));
t24 = sin(qJ(1));
t37 = t20 * t24;
t36 = t20 * t26;
t35 = cos(pkin(6));
t34 = r_i_i_C(3) + pkin(10) - pkin(9);
t33 = t20 * t38;
t23 = sin(qJ(2));
t27 = cos(qJ(2));
t31 = t35 * t38;
t12 = t23 * t31 + t24 * t27;
t4 = t12 * t26 - t22 * t33;
t32 = t24 * t35;
t3 = t12 * t22 + t26 * t33;
t14 = -t23 * t32 + t38 * t27;
t13 = t38 * t23 + t27 * t32;
t11 = t24 * t23 - t27 * t31;
t10 = t35 * t22 + t23 * t36;
t9 = t20 * t23 * t22 - t35 * t26;
t8 = t14 * t26 + t22 * t37;
t7 = t14 * t22 - t24 * t36;
t2 = t7 * t21 + t8 * t25;
t1 = -t8 * t21 + t7 * t25;
t5 = [-t24 * pkin(1) - t12 * pkin(2) + pkin(8) * t33 + t34 * t11 - t29 * t4 - t30 * t3, -t13 * t40 - t34 * t14, -t29 * t7 + t30 * t8, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t38 * pkin(1) + t14 * pkin(2) + pkin(8) * t37 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t7 * qJ(4) - t34 * t13 + t39 * t8, -t11 * t40 - t34 * t12, -t29 * t3 + t30 * t4, t3 (-t4 * t21 + t3 * t25) * r_i_i_C(1) + (-t3 * t21 - t4 * t25) * r_i_i_C(2), 0; 0 (-t34 * t23 + t40 * t27) * t20, t30 * t10 - t29 * t9, t9 (-t10 * t21 + t9 * t25) * r_i_i_C(1) + (-t10 * t25 - t9 * t21) * r_i_i_C(2), 0;];
Ja_transl  = t5;
