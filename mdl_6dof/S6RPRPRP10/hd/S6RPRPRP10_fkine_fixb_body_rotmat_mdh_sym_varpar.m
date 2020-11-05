% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRP10 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:38
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRPRP10_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP10_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:38:49
	% EndTime: 2020-11-04 21:38:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:38:49
	% EndTime: 2020-11-04 21:38:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t48 = cos(qJ(1));
	t47 = sin(qJ(1));
	t1 = [t48, -t47, 0, 0; t47, t48, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:38:49
	% EndTime: 2020-11-04 21:38:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t50 = cos(qJ(1));
	t49 = sin(qJ(1));
	t1 = [0, -t50, t49, t50 * pkin(1) + t49 * qJ(2) + 0; 0, -t49, -t50, t49 * pkin(1) - t50 * qJ(2) + 0; 1, 0, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:38:49
	% EndTime: 2020-11-04 21:38:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t55 = pkin(1) + pkin(7);
	t54 = cos(qJ(1));
	t53 = cos(qJ(3));
	t52 = sin(qJ(1));
	t51 = sin(qJ(3));
	t1 = [t52 * t51, t52 * t53, t54, t52 * qJ(2) + t55 * t54 + 0; -t54 * t51, -t54 * t53, t52, -t54 * qJ(2) + t55 * t52 + 0; t53, -t51, 0, pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:38:49
	% EndTime: 2020-11-04 21:38:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (17->14), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->7)
	t61 = pkin(1) + pkin(7);
	t60 = cos(qJ(1));
	t59 = cos(qJ(3));
	t58 = sin(qJ(1));
	t57 = sin(qJ(3));
	t56 = -t57 * pkin(3) + t59 * qJ(4) - qJ(2);
	t1 = [t60, -t58 * t57, -t58 * t59, -t56 * t58 + t61 * t60 + 0; t58, t60 * t57, t60 * t59, t56 * t60 + t61 * t58 + 0; 0, -t59, t57, t59 * pkin(3) + t57 * qJ(4) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:38:49
	% EndTime: 2020-11-04 21:38:49
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (24->17), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->13)
	t65 = sin(qJ(1));
	t67 = cos(qJ(3));
	t73 = t65 * t67;
	t63 = sin(qJ(5));
	t68 = cos(qJ(1));
	t72 = t68 * t63;
	t66 = cos(qJ(5));
	t71 = t68 * t66;
	t64 = sin(qJ(3));
	t69 = pkin(3) + pkin(8);
	t70 = t67 * qJ(4) - t69 * t64 - qJ(2);
	t62 = pkin(1) + pkin(4) + pkin(7);
	t1 = [-t63 * t73 + t71, -t66 * t73 - t72, t65 * t64, t62 * t68 - t70 * t65 + 0; t65 * t66 + t67 * t72, -t65 * t63 + t67 * t71, -t68 * t64, t62 * t65 + t70 * t68 + 0; t64 * t63, t64 * t66, t67, t64 * qJ(4) + t69 * t67 + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:38:49
	% EndTime: 2020-11-04 21:38:49
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (35->23), mult. (36->25), div. (0->0), fcn. (49->6), ass. (0->14)
	t79 = sin(qJ(1));
	t81 = cos(qJ(3));
	t87 = t79 * t81;
	t77 = sin(qJ(5));
	t82 = cos(qJ(1));
	t86 = t82 * t77;
	t80 = cos(qJ(5));
	t85 = t82 * t80;
	t75 = -pkin(5) * t77 + qJ(6) * t80 - qJ(4);
	t78 = sin(qJ(3));
	t83 = pkin(3) + pkin(8);
	t84 = t75 * t81 + t78 * t83 + qJ(2);
	t74 = pkin(5) * t80 + qJ(6) * t77 + pkin(1) + pkin(4) + pkin(7);
	t1 = [-t77 * t87 + t85, t79 * t78, t80 * t87 + t86, t74 * t82 + t79 * t84 + 0; t79 * t80 + t81 * t86, -t82 * t78, t77 * t79 - t81 * t85, t74 * t79 - t82 * t84 + 0; t78 * t77, t81, -t78 * t80, -t75 * t78 + t81 * t83 + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end