% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRP2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:05
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PRRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:05:36
	% EndTime: 2020-11-04 20:05:36
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:05:36
	% EndTime: 2020-11-04 20:05:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t39 = cos(pkin(8));
	t38 = sin(pkin(8));
	t1 = [t39, -t38, 0, 0; t38, t39, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:05:36
	% EndTime: 2020-11-04 20:05:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t42 = pkin(8) + qJ(2);
	t41 = cos(t42);
	t40 = sin(t42);
	t1 = [t41, -t40, 0, cos(pkin(8)) * pkin(1) + 0; t40, t41, 0, sin(pkin(8)) * pkin(1) + 0; 0, 0, 1, pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:05:36
	% EndTime: 2020-11-04 20:05:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (18->10), mult. (4->4), div. (0->0), fcn. (8->6), ass. (0->5)
	t46 = pkin(8) + qJ(2);
	t45 = qJ(3) + t46;
	t44 = cos(t45);
	t43 = sin(t45);
	t1 = [t44, -t43, 0, pkin(2) * cos(t46) + cos(pkin(8)) * pkin(1) + 0; t43, t44, 0, pkin(2) * sin(t46) + sin(pkin(8)) * pkin(1) + 0; 0, 0, 1, pkin(6) + pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:05:36
	% EndTime: 2020-11-04 20:05:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (36->16), mult. (12->12), div. (0->0), fcn. (20->8), ass. (0->7)
	t50 = pkin(8) + qJ(2);
	t52 = cos(qJ(4));
	t51 = sin(qJ(4));
	t49 = qJ(3) + t50;
	t48 = cos(t49);
	t47 = sin(t49);
	t1 = [t48 * t52, -t48 * t51, t47, t48 * pkin(3) + t47 * pkin(7) + pkin(2) * cos(t50) + cos(pkin(8)) * pkin(1) + 0; t47 * t52, -t47 * t51, -t48, t47 * pkin(3) - t48 * pkin(7) + pkin(2) * sin(t50) + sin(pkin(8)) * pkin(1) + 0; t51, t52, 0, pkin(6) + pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:05:36
	% EndTime: 2020-11-04 20:05:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (41->19), mult. (18->16), div. (0->0), fcn. (26->8), ass. (0->8)
	t57 = pkin(8) + qJ(2);
	t59 = cos(qJ(4));
	t58 = sin(qJ(4));
	t56 = qJ(3) + t57;
	t55 = cos(t56);
	t54 = sin(t56);
	t53 = pkin(4) * t59 + qJ(5) * t58 + pkin(3);
	t1 = [t55 * t59, t54, t55 * t58, t53 * t55 + cos(pkin(8)) * pkin(1) + pkin(2) * cos(t57) + t54 * pkin(7) + 0; t54 * t59, -t55, t54 * t58, t53 * t54 + sin(pkin(8)) * pkin(1) + pkin(2) * sin(t57) - t55 * pkin(7) + 0; t58, 0, -t59, t58 * pkin(4) - t59 * qJ(5) + pkin(5) + pkin(6) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end