% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRPR6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:26
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPPRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:26:27
	% EndTime: 2020-11-04 21:26:27
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:26:27
	% EndTime: 2020-11-04 21:26:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t39 = cos(qJ(1));
	t38 = sin(qJ(1));
	t1 = [t39, -t38, 0, 0; t38, t39, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:26:27
	% EndTime: 2020-11-04 21:26:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t41 = cos(qJ(1));
	t40 = sin(qJ(1));
	t1 = [0, -t41, t40, t41 * pkin(1) + t40 * qJ(2) + 0; 0, -t40, -t41, t40 * pkin(1) - t41 * qJ(2) + 0; 1, 0, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:26:27
	% EndTime: 2020-11-04 21:26:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->4)
	t44 = cos(qJ(1));
	t43 = sin(qJ(1));
	t42 = pkin(1) + qJ(3);
	t1 = [0, t43, t44, t43 * qJ(2) + t42 * t44 + 0; 0, -t44, t43, -t44 * qJ(2) + t42 * t43 + 0; 1, 0, 0, pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:26:27
	% EndTime: 2020-11-04 21:26:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (13->11), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->7)
	t50 = cos(qJ(1));
	t49 = cos(qJ(4));
	t48 = sin(qJ(1));
	t47 = sin(qJ(4));
	t46 = pkin(1) + qJ(3);
	t45 = pkin(7) - qJ(2);
	t1 = [t50 * t47, t50 * t49, -t48, -t45 * t48 + t46 * t50 + 0; t48 * t47, t48 * t49, t50, t45 * t50 + t46 * t48 + 0; t49, -t47, 0, pkin(3) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:26:27
	% EndTime: 2020-11-04 21:26:27
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (23->19), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->7)
	t56 = cos(qJ(1));
	t55 = cos(qJ(4));
	t54 = sin(qJ(1));
	t53 = sin(qJ(4));
	t52 = pkin(7) - qJ(2);
	t51 = t53 * pkin(4) - t55 * qJ(5) + pkin(1) + qJ(3);
	t1 = [-t54, -t56 * t53, -t56 * t55, t51 * t56 - t52 * t54 + 0; t56, -t54 * t53, -t54 * t55, t51 * t54 + t52 * t56 + 0; 0, -t55, t53, t55 * pkin(4) + t53 * qJ(5) + pkin(2) + pkin(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:26:27
	% EndTime: 2020-11-04 21:26:27
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (26->18), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->13)
	t61 = sin(qJ(1));
	t63 = cos(qJ(4));
	t68 = t61 * t63;
	t59 = sin(qJ(6));
	t64 = cos(qJ(1));
	t67 = t64 * t59;
	t62 = cos(qJ(6));
	t66 = t64 * t62;
	t65 = pkin(4) + pkin(8);
	t60 = sin(qJ(4));
	t58 = pkin(5) + pkin(7) - qJ(2);
	t57 = -t63 * qJ(5) + t65 * t60 + pkin(1) + qJ(3);
	t1 = [-t61 * t62 - t63 * t67, t61 * t59 - t63 * t66, t64 * t60, t57 * t64 - t58 * t61 + 0; -t59 * t68 + t66, -t62 * t68 - t67, t61 * t60, t57 * t61 + t58 * t64 + 0; t60 * t59, t60 * t62, t63, t60 * qJ(5) + t65 * t63 + pkin(2) + pkin(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end