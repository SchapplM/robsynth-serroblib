% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR5_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR5_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:18:17
% EndTime: 2019-02-26 22:18:17
% DurationCPUTime: 0.15s
% Computational Cost: add. (122->28), mult. (139->40), div. (0->0), fcn. (150->8), ass. (0->27)
t22 = sin(qJ(5));
t25 = cos(qJ(5));
t48 = r_i_i_C(1) * t22 + r_i_i_C(2) * t25;
t21 = qJ(2) + qJ(3);
t18 = sin(t21);
t19 = cos(t21);
t36 = pkin(3) + pkin(9) + r_i_i_C(3);
t47 = t18 * qJ(4) + t36 * t19;
t45 = (qJ(4) + t48) * t19;
t20 = cos(qJ(2)) * pkin(2);
t44 = t20 + pkin(1) + t47;
t41 = pkin(4) + pkin(8) + pkin(7);
t24 = sin(qJ(1));
t40 = t24 * t18;
t39 = t24 * t25;
t26 = cos(qJ(1));
t38 = t26 * t18;
t35 = t45 * t24;
t34 = t45 * t26;
t30 = t36 * t18;
t29 = t48 * t18 + t47;
t28 = -sin(qJ(2)) * pkin(2) - t30;
t4 = -t22 * t40 + t25 * t26;
t3 = t18 * t39 + t22 * t26;
t2 = t22 * t38 + t39;
t1 = -t22 * t24 + t25 * t38;
t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t44 * t24 + t41 * t26, t28 * t26 + t34, -t26 * t30 + t34, t38, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t41 * t24 + t44 * t26, t28 * t24 + t35, -t24 * t30 + t35, t40, r_i_i_C(1) * t3 + r_i_i_C(2) * t4, 0; 0, t20 + t29, t29, -t19 (-r_i_i_C(1) * t25 + r_i_i_C(2) * t22) * t19, 0;];
Ja_transl  = t5;
