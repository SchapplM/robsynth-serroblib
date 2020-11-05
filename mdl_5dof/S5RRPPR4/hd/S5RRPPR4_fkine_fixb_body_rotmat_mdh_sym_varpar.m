% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPPR4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:30
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:30:19
	% EndTime: 2020-11-04 20:30:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:30:19
	% EndTime: 2020-11-04 20:30:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t42 = cos(qJ(1));
	t41 = sin(qJ(1));
	t1 = [t42, -t41, 0, 0; t41, t42, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:30:19
	% EndTime: 2020-11-04 20:30:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t45 = qJ(1) + qJ(2);
	t44 = cos(t45);
	t43 = sin(t45);
	t1 = [t44, -t43, 0, cos(qJ(1)) * pkin(1) + 0; t43, t44, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:30:19
	% EndTime: 2020-11-04 20:30:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (17->10), mult. (6->6), div. (0->0), fcn. (10->4), ass. (0->4)
	t48 = qJ(1) + qJ(2);
	t47 = cos(t48);
	t46 = sin(t48);
	t1 = [t47, 0, t46, t47 * pkin(2) + t46 * qJ(3) + cos(qJ(1)) * pkin(1) + 0; t46, 0, -t47, t46 * pkin(2) - t47 * qJ(3) + sin(qJ(1)) * pkin(1) + 0; 0, 1, 0, pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:30:19
	% EndTime: 2020-11-04 20:30:19
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (27->14), mult. (14->10), div. (0->0), fcn. (22->6), ass. (0->9)
	t56 = pkin(2) + pkin(3);
	t55 = cos(pkin(8));
	t54 = sin(pkin(8));
	t53 = qJ(1) + qJ(2);
	t52 = cos(t53);
	t51 = sin(t53);
	t50 = t51 * t55 - t52 * t54;
	t49 = -t51 * t54 - t52 * t55;
	t1 = [-t49, t50, 0, t56 * t52 + t51 * qJ(3) + cos(qJ(1)) * pkin(1) + 0; t50, t49, 0, t56 * t51 - t52 * qJ(3) + sin(qJ(1)) * pkin(1) + 0; 0, 0, -1, -qJ(4) + pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:30:19
	% EndTime: 2020-11-04 20:30:19
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (49->22), mult. (26->18), div. (0->0), fcn. (40->10), ass. (0->14)
	t70 = cos(pkin(8));
	t65 = qJ(1) + qJ(2);
	t69 = pkin(2) + pkin(3);
	t68 = cos(qJ(5));
	t67 = sin(qJ(5));
	t66 = sin(pkin(8));
	t64 = cos(t65);
	t63 = sin(t65);
	t62 = pkin(8) - t65;
	t61 = cos(t62);
	t60 = sin(t62);
	t58 = t63 * t66 + t64 * t70;
	t57 = -t63 * t70 + t64 * t66;
	t1 = [t58 * t68, -t58 * t67, t57, pkin(4) * t61 + pkin(7) * t60 + t69 * t64 + t63 * qJ(3) + cos(qJ(1)) * pkin(1) + 0; -t57 * t68, t57 * t67, t58, pkin(7) * t61 - pkin(4) * t60 + t69 * t63 - t64 * qJ(3) + sin(qJ(1)) * pkin(1) + 0; -t67, -t68, 0, -qJ(4) + pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end