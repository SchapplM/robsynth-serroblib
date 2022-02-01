% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPR2
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRPR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
JaD_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:31:18
	% EndTime: 2022-01-20 11:31:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:31:18
	% EndTime: 2022-01-20 11:31:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:31:18
	% EndTime: 2022-01-20 11:31:18
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (22->6), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->8)
	t43 = pkin(1) * qJD(1);
	t40 = qJ(1) + qJ(2);
	t37 = sin(t40);
	t38 = cos(t40);
	t39 = qJD(1) + qJD(2);
	t42 = (-r_i_i_C(1) * t38 + r_i_i_C(2) * t37) * t39;
	t41 = (-r_i_i_C(1) * t37 - r_i_i_C(2) * t38) * t39;
	t1 = [-cos(qJ(1)) * t43 + t42, t42, 0, 0, 0; -sin(qJ(1)) * t43 + t41, t41, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:31:18
	% EndTime: 2022-01-20 11:31:18
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (68->10), mult. (36->12), div. (0->0), fcn. (18->6), ass. (0->13)
	t48 = qJD(1) + qJD(2);
	t55 = pkin(2) * t48;
	t54 = pkin(1) * qJD(1);
	t49 = qJ(1) + qJ(2);
	t47 = qJ(3) + t49;
	t42 = sin(t47);
	t43 = cos(t47);
	t44 = qJD(3) + t48;
	t53 = (-r_i_i_C(1) * t43 + r_i_i_C(2) * t42) * t44;
	t52 = (-r_i_i_C(1) * t42 - r_i_i_C(2) * t43) * t44;
	t51 = -cos(t49) * t55 + t53;
	t50 = -sin(t49) * t55 + t52;
	t1 = [-cos(qJ(1)) * t54 + t51, t51, t53, 0, 0; -sin(qJ(1)) * t54 + t50, t50, t52, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:31:18
	% EndTime: 2022-01-20 11:31:18
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (110->13), mult. (48->14), div. (0->0), fcn. (24->8), ass. (0->14)
	t54 = qJD(1) + qJD(2);
	t61 = pkin(2) * t54;
	t60 = pkin(1) * qJD(1);
	t55 = qJ(1) + qJ(2);
	t53 = qJ(3) + t55;
	t47 = pkin(9) + t53;
	t45 = sin(t47);
	t46 = cos(t47);
	t50 = qJD(3) + t54;
	t59 = (-pkin(3) * cos(t53) - r_i_i_C(1) * t46 + t45 * r_i_i_C(2)) * t50;
	t58 = (-pkin(3) * sin(t53) - r_i_i_C(1) * t45 - r_i_i_C(2) * t46) * t50;
	t57 = -cos(t55) * t61 + t59;
	t56 = -sin(t55) * t61 + t58;
	t1 = [-cos(qJ(1)) * t60 + t57, t57, t59, 0, 0; -sin(qJ(1)) * t60 + t56, t56, t58, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:31:18
	% EndTime: 2022-01-20 11:31:18
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (305->26), mult. (160->38), div. (0->0), fcn. (96->10), ass. (0->26)
	t78 = r_i_i_C(3) + pkin(8);
	t57 = qJ(1) + qJ(2);
	t55 = qJ(3) + t57;
	t49 = pkin(9) + t55;
	t47 = sin(t49);
	t48 = cos(t49);
	t59 = cos(qJ(5));
	t70 = qJD(5) * t59;
	t56 = qJD(1) + qJD(2);
	t52 = qJD(3) + t56;
	t58 = sin(qJ(5));
	t74 = t52 * t58;
	t77 = t47 * t70 + t48 * t74;
	t76 = pkin(2) * t56;
	t73 = t52 * t59;
	t72 = pkin(1) * qJD(1);
	t71 = qJD(5) * t58;
	t69 = t47 * t74;
	t67 = t47 * t71;
	t65 = -r_i_i_C(1) * t59 - pkin(4);
	t64 = (-r_i_i_C(1) * t58 - r_i_i_C(2) * t59) * qJD(5);
	t63 = r_i_i_C(1) * t67 + (-pkin(3) * cos(t55) + t65 * t48 - t78 * t47) * t52 + t77 * r_i_i_C(2);
	t62 = -cos(t57) * t76 + t63;
	t61 = r_i_i_C(2) * t69 + (-pkin(3) * sin(t55) + t65 * t47) * t52 + (t78 * t52 + t64) * t48;
	t60 = -sin(t57) * t76 + t61;
	t1 = [-cos(qJ(1)) * t72 + t62, t62, t63, 0, (t47 * t73 + t48 * t71) * r_i_i_C(2) + (-t48 * t70 + t69) * r_i_i_C(1); -sin(qJ(1)) * t72 + t60, t60, t61, 0, (-t48 * t73 + t67) * r_i_i_C(2) - t77 * r_i_i_C(1); 0, 0, 0, 0, t64;];
	JaD_transl = t1;
end