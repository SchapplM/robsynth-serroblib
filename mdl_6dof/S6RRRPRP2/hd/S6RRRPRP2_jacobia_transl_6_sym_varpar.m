% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRP2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:10:09
% EndTime: 2019-02-26 22:10:09
% DurationCPUTime: 0.15s
% Computational Cost: add. (247->33), mult. (206->41), div. (0->0), fcn. (229->10), ass. (0->31)
t28 = sin(qJ(5));
t30 = cos(qJ(5));
t38 = r_i_i_C(3) + qJ(6);
t50 = pkin(5) + r_i_i_C(1);
t57 = t38 * t28 + t50 * t30;
t56 = pkin(9) + r_i_i_C(2);
t27 = qJ(2) + qJ(3);
t23 = pkin(10) + t27;
t19 = cos(t23);
t54 = t19 * t56;
t18 = sin(t23);
t53 = t56 * t18 + t19 * pkin(4) + pkin(3) * cos(t27);
t32 = (-pkin(4) - t57) * t18 - pkin(3) * sin(t27);
t25 = cos(qJ(2)) * pkin(2);
t51 = pkin(1) + t25 + t53;
t29 = sin(qJ(1));
t48 = t29 * t54;
t42 = t29 * t28;
t41 = t29 * t30;
t31 = cos(qJ(1));
t40 = t31 * t28;
t39 = t31 * t30;
t37 = t31 * t54;
t35 = t57 * t19 + t53;
t33 = -sin(qJ(2)) * pkin(2) + t32;
t26 = -qJ(4) - pkin(8) - pkin(7);
t4 = t19 * t39 + t42;
t3 = t19 * t40 - t41;
t2 = t19 * t41 - t40;
t1 = t19 * t42 + t39;
t5 = [-t38 * t1 - t50 * t2 - t31 * t26 - t29 * t51, t33 * t31 + t37, t32 * t31 + t37, t29, -t50 * t3 + t38 * t4, t3; -t29 * t26 + t38 * t3 + t31 * t51 + t50 * t4, t33 * t29 + t48, t32 * t29 + t48, -t31, -t50 * t1 + t38 * t2, t1; 0, t25 + t35, t35, 0 (-t50 * t28 + t38 * t30) * t18, t18 * t28;];
Ja_transl  = t5;
