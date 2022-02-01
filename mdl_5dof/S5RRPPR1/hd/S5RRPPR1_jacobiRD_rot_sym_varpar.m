% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: S5RRPPR1_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRPPR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
JRD_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:52:15
	% EndTime: 2022-01-20 09:52:15
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:52:15
	% EndTime: 2022-01-20 09:52:15
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:52:15
	% EndTime: 2022-01-20 09:52:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->8), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t47 = qJD(1) + qJD(2);
	t48 = qJ(1) + qJ(2);
	t49 = t47 * cos(t48);
	t44 = t47 * sin(t48);
	t1 = [-t49, -t49, 0, 0, 0; -t44, -t44, 0, 0, 0; 0, 0, 0, 0, 0; t44, t44, 0, 0, 0; -t49, -t49, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:52:15
	% EndTime: 2022-01-20 09:52:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (30->8), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t53 = qJ(1) + qJ(2) + pkin(8);
	t54 = qJD(1) + qJD(2);
	t55 = t54 * cos(t53);
	t50 = t54 * sin(t53);
	t1 = [-t55, -t55, 0, 0, 0; -t50, -t50, 0, 0, 0; 0, 0, 0, 0, 0; t50, t50, 0, 0, 0; -t55, -t55, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:52:15
	% EndTime: 2022-01-20 09:52:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (42->8), mult. (20->8), div. (0->0), fcn. (20->4), ass. (0->13)
	t59 = qJ(1) + qJ(2) + pkin(8);
	t57 = sin(t59);
	t60 = qJD(1) + qJD(2);
	t67 = t60 * t57;
	t66 = t60 * sin(pkin(9));
	t65 = t60 * cos(pkin(9));
	t64 = t57 * t65;
	t58 = cos(t59);
	t63 = t58 * t65;
	t56 = t60 * t58;
	t55 = t58 * t66;
	t54 = t57 * t66;
	t1 = [-t63, -t63, 0, 0, 0; -t64, -t64, 0, 0, 0; 0, 0, 0, 0, 0; t55, t55, 0, 0, 0; t54, t54, 0, 0, 0; 0, 0, 0, 0, 0; -t67, -t67, 0, 0, 0; t56, t56, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:52:15
	% EndTime: 2022-01-20 09:52:15
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (114->14), mult. (54->12), div. (0->0), fcn. (54->4), ass. (0->16)
	t81 = qJ(1) + qJ(2) + pkin(8);
	t77 = sin(t81);
	t83 = qJD(1) + qJD(2);
	t86 = t83 * t77;
	t78 = cos(t81);
	t76 = t83 * t78;
	t82 = pkin(9) + qJ(5);
	t79 = sin(t82);
	t85 = qJD(5) * t79;
	t80 = cos(t82);
	t84 = qJD(5) * t80;
	t75 = -t80 * t76 + t77 * t85;
	t74 = t79 * t76 + t77 * t84;
	t73 = t78 * t85 + t80 * t86;
	t72 = -t78 * t84 + t79 * t86;
	t1 = [t75, t75, 0, 0, t72; -t73, -t73, 0, 0, -t74; 0, 0, 0, 0, -t85; t74, t74, 0, 0, t73; t72, t72, 0, 0, t75; 0, 0, 0, 0, -t84; -t86, -t86, 0, 0, 0; t76, t76, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
end