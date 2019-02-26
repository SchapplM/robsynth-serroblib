% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRP4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:11:05
% EndTime: 2019-02-26 22:11:05
% DurationCPUTime: 0.14s
% Computational Cost: add. (179->33), mult. (214->45), div. (0->0), fcn. (238->8), ass. (0->31)
t22 = qJ(2) + qJ(3);
t19 = sin(t22);
t20 = cos(t22);
t38 = pkin(3) + pkin(9) + r_i_i_C(2);
t51 = t19 * qJ(4) + t38 * t20;
t47 = pkin(5) + r_i_i_C(1);
t49 = t20 * t47;
t21 = cos(qJ(2)) * pkin(2);
t48 = t21 + pkin(1) + t51;
t46 = pkin(4) + pkin(8) + pkin(7);
t23 = sin(qJ(5));
t25 = sin(qJ(1));
t44 = t25 * t23;
t26 = cos(qJ(5));
t43 = t25 * t26;
t27 = cos(qJ(1));
t42 = t27 * t23;
t41 = t27 * t26;
t40 = r_i_i_C(3) + qJ(6);
t39 = qJ(4) * t20;
t37 = t25 * t39 + t44 * t49;
t36 = t27 * t39 + t42 * t49;
t33 = t40 * t26;
t31 = -t38 * t19 - t20 * t33;
t30 = (t23 * t47 - t33) * t19 + t51;
t29 = -sin(qJ(2)) * pkin(2) + t31;
t4 = -t19 * t44 + t41;
t3 = t19 * t43 + t42;
t2 = t19 * t42 + t43;
t1 = -t19 * t41 + t44;
t5 = [-t48 * t25 + t46 * t27 + t40 * t3 + t47 * t4, t29 * t27 + t36, t31 * t27 + t36, t27 * t19, -t47 * t1 + t40 * t2, t1; t40 * t1 + t47 * t2 + t46 * t25 + t48 * t27, t29 * t25 + t37, t31 * t25 + t37, t25 * t19, t47 * t3 - t40 * t4, -t3; 0, t21 + t30, t30, -t20 (-t40 * t23 - t47 * t26) * t20, t20 * t26;];
Ja_transl  = t5;
