% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP12
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

function Ja_transl = S6RRRPRP12_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP12_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:15:31
% EndTime: 2019-02-26 22:15:31
% DurationCPUTime: 0.21s
% Computational Cost: add. (275->61), mult. (691->99), div. (0->0), fcn. (888->10), ass. (0->40)
t35 = sin(qJ(2));
t36 = sin(qJ(1));
t39 = cos(qJ(2));
t48 = cos(pkin(6));
t58 = cos(qJ(1));
t43 = t48 * t58;
t26 = t35 * t43 + t36 * t39;
t34 = sin(qJ(3));
t38 = cos(qJ(3));
t32 = sin(pkin(6));
t46 = t32 * t58;
t15 = t26 * t34 + t38 * t46;
t25 = t36 * t35 - t39 * t43;
t33 = sin(qJ(5));
t37 = cos(qJ(5));
t63 = t15 * t33 + t25 * t37;
t62 = t15 * t37 - t25 * t33;
t47 = pkin(3) + pkin(10) + r_i_i_C(2);
t61 = qJ(4) * t34 + t47 * t38 + pkin(2);
t60 = pkin(4) + pkin(9);
t59 = pkin(5) + r_i_i_C(1);
t55 = t32 * t36;
t54 = t32 * t38;
t53 = t33 * t34;
t52 = t33 * t39;
t51 = t34 * t37;
t50 = t37 * t39;
t49 = r_i_i_C(3) + qJ(6);
t44 = t36 * t48;
t42 = t26 * t38 - t34 * t46;
t40 = t59 * t33 - t49 * t37 + qJ(4);
t28 = -t35 * t44 + t58 * t39;
t27 = t58 * t35 + t39 * t44;
t23 = t32 * t35 * t34 - t48 * t38;
t20 = t28 * t38 + t34 * t55;
t19 = t28 * t34 - t36 * t54;
t13 = t23 * t37 + t32 * t52;
t6 = t19 * t33 + t27 * t37;
t5 = -t19 * t37 + t27 * t33;
t1 = [-t36 * pkin(1) - t26 * pkin(2) + pkin(8) * t46 - t15 * qJ(4) - t60 * t25 - t47 * t42 + t49 * t62 - t59 * t63, t49 * (t27 * t51 + t28 * t33) + t60 * t28 + t59 * (-t27 * t53 + t28 * t37) - t61 * t27, -t47 * t19 + t40 * t20, t19, t49 * t6 - t59 * t5, t5; t58 * pkin(1) + t28 * pkin(2) + pkin(8) * t55 + t19 * qJ(4) + t47 * t20 + t60 * t27 + t49 * t5 + t59 * t6, t59 * (-t25 * t53 + t26 * t37) + t49 * (t25 * t51 + t26 * t33) + t60 * t26 - t61 * t25, -t47 * t15 + t40 * t42, t15, t49 * t63 + t59 * t62, -t62; 0 (t59 * (t34 * t52 + t35 * t37) + t49 * (t33 * t35 - t34 * t50) + t60 * t35 + t61 * t39) * t32, -t47 * t23 + t40 * (t48 * t34 + t35 * t54) t23, -t49 * (-t23 * t33 + t32 * t50) + t59 * t13, -t13;];
Ja_transl  = t1;
