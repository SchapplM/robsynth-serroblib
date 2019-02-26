% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Ja_transl = S6RRRPRP2_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP2_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:10:04
% EndTime: 2019-02-26 22:10:04
% DurationCPUTime: 0.15s
% Computational Cost: add. (165->32), mult. (131->41), div. (0->0), fcn. (141->10), ass. (0->30)
t24 = qJ(2) + qJ(3);
t20 = pkin(10) + t24;
t17 = sin(t20);
t18 = cos(t20);
t25 = sin(qJ(5));
t43 = r_i_i_C(2) * t25;
t51 = pkin(9) + r_i_i_C(3);
t52 = t17 * t43 + t18 * t51;
t49 = t51 * t17 + t18 * pkin(4) + pkin(3) * cos(t24);
t27 = cos(qJ(5));
t44 = r_i_i_C(1) * t27;
t30 = (-pkin(4) - t44) * t17 - pkin(3) * sin(t24);
t22 = cos(qJ(2)) * pkin(2);
t47 = pkin(1) + t22 + t49;
t28 = cos(qJ(1));
t40 = t25 * t28;
t26 = sin(qJ(1));
t39 = t26 * t25;
t38 = t26 * t27;
t37 = t27 * t28;
t36 = t52 * t26;
t34 = t52 * t28;
t31 = -sin(qJ(2)) * pkin(2) + t30;
t29 = (-t43 + t44) * t18 + t49;
t23 = -qJ(4) - pkin(8) - pkin(7);
t4 = t18 * t37 + t39;
t3 = -t18 * t40 + t38;
t2 = -t18 * t38 + t40;
t1 = t18 * t39 + t37;
t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t23 * t28 - t47 * t26, t31 * t28 + t34, t30 * t28 + t34, t26, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t26 * t23 + t47 * t28, t31 * t26 + t36, t30 * t26 + t36, -t28, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0; 0, t22 + t29, t29, 0 (-r_i_i_C(1) * t25 - r_i_i_C(2) * t27) * t17, 0;];
Ja_transl  = t5;
