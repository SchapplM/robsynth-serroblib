% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPP9_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP9_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:30:05
% EndTime: 2019-02-26 22:30:06
% DurationCPUTime: 0.20s
% Computational Cost: add. (242->58), mult. (613->99), div. (0->0), fcn. (787->10), ass. (0->39)
t34 = sin(qJ(2));
t35 = sin(qJ(1));
t38 = cos(qJ(2));
t39 = cos(qJ(1));
t46 = cos(pkin(6));
t43 = t39 * t46;
t25 = t34 * t43 + t35 * t38;
t33 = sin(qJ(3));
t37 = cos(qJ(3));
t31 = sin(pkin(6));
t51 = t31 * t39;
t14 = t25 * t37 - t33 * t51;
t24 = t35 * t34 - t38 * t43;
t32 = sin(qJ(4));
t36 = cos(qJ(4));
t1 = t14 * t32 - t24 * t36;
t60 = t14 * t36 + t24 * t32;
t57 = pkin(10) + r_i_i_C(1);
t59 = pkin(3) * t37 + t57 * t33 + pkin(2);
t47 = r_i_i_C(3) + qJ(5);
t58 = pkin(4) - r_i_i_C(2);
t40 = t47 * t32 + t58 * t36 + pkin(3);
t54 = t31 * t35;
t53 = t31 * t37;
t52 = t31 * t38;
t50 = t32 * t37;
t49 = t36 * t37;
t48 = t37 * t38;
t44 = t35 * t46;
t42 = -t25 * t33 - t37 * t51;
t27 = -t34 * t44 + t39 * t38;
t26 = t39 * t34 + t38 * t44;
t23 = t46 * t33 + t34 * t53;
t18 = t27 * t37 + t33 * t54;
t17 = t27 * t33 - t35 * t53;
t11 = t23 * t32 + t36 * t52;
t6 = t18 * t36 + t26 * t32;
t5 = t18 * t32 - t26 * t36;
t2 = [-t35 * pkin(1) - t25 * pkin(2) - t14 * pkin(3) + pkin(8) * t51 - t24 * pkin(9) - t47 * t1 + t57 * t42 - t58 * t60, t27 * pkin(9) + t47 * (-t26 * t50 - t27 * t36) + t58 * (-t26 * t49 + t27 * t32) - t59 * t26, -t40 * t17 + t57 * t18, t47 * t6 - t58 * t5, t5, 0; t39 * pkin(1) + t27 * pkin(2) + t18 * pkin(3) + pkin(8) * t54 + t26 * pkin(9) + t57 * t17 + t47 * t5 + t58 * t6, t25 * pkin(9) + t58 * (-t24 * t49 + t25 * t32) + t47 * (-t24 * t50 - t25 * t36) - t59 * t24, t57 * t14 + t40 * t42, -t58 * t1 + t47 * t60, t1, 0; 0 (t58 * (t32 * t34 + t36 * t48) + t47 * (t32 * t48 - t34 * t36) + pkin(9) * t34 + t59 * t38) * t31, t57 * t23 + t40 * (-t31 * t34 * t33 + t46 * t37) t47 * (t23 * t36 - t32 * t52) - t58 * t11, t11, 0;];
Ja_transl  = t2;
