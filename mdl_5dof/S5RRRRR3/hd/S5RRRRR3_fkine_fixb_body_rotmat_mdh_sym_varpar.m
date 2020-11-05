% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRR3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:48
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:48:53
	% EndTime: 2020-11-04 20:48:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:48:53
	% EndTime: 2020-11-04 20:48:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t47 = cos(qJ(1));
	t46 = sin(qJ(1));
	t1 = [t47, -t46, 0, 0; t46, t47, 0, 0; 0, 0, 1, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:48:53
	% EndTime: 2020-11-04 20:48:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (4->4), mult. (4->4), div. (0->0), fcn. (12->4), ass. (0->5)
	t51 = cos(qJ(1));
	t50 = cos(qJ(2));
	t49 = sin(qJ(1));
	t48 = sin(qJ(2));
	t1 = [t51 * t50, -t51 * t48, t49, 0; t49 * t50, -t49 * t48, -t51, 0; t48, t50, 0, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:48:53
	% EndTime: 2020-11-04 20:48:53
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (13->8), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->7)
	t58 = pkin(1) * cos(qJ(2));
	t57 = cos(qJ(1));
	t55 = sin(qJ(1));
	t54 = qJ(2) + qJ(3);
	t53 = cos(t54);
	t52 = sin(t54);
	t1 = [t57 * t53, -t57 * t52, t55, t57 * t58 + 0; t55 * t53, -t55 * t52, -t57, t55 * t58 + 0; t52, t53, 0, sin(qJ(2)) * pkin(1) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:48:53
	% EndTime: 2020-11-04 20:48:53
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (31->15), mult. (31->20), div. (0->0), fcn. (44->8), ass. (0->13)
	t62 = sin(qJ(4));
	t63 = sin(qJ(1));
	t71 = t63 * t62;
	t64 = cos(qJ(4));
	t70 = t63 * t64;
	t66 = cos(qJ(1));
	t69 = t66 * t62;
	t68 = t66 * t64;
	t61 = qJ(2) + qJ(3);
	t59 = sin(t61);
	t60 = cos(t61);
	t67 = pkin(1) * cos(qJ(2)) + pkin(2) * t60 + pkin(5) * t59;
	t1 = [t60 * t68 + t71, -t60 * t69 + t70, t66 * t59, t67 * t66 + 0; t60 * t70 - t69, -t60 * t71 - t68, t63 * t59, t67 * t63 + 0; t59 * t64, -t59 * t62, -t60, t59 * pkin(2) - t60 * pkin(5) + sin(qJ(2)) * pkin(1) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:48:53
	% EndTime: 2020-11-04 20:48:53
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (46->19), mult. (38->24), div. (0->0), fcn. (51->10), ass. (0->16)
	t88 = pkin(3) * sin(qJ(4));
	t77 = qJ(4) + qJ(5);
	t73 = sin(t77);
	t80 = sin(qJ(1));
	t87 = t80 * t73;
	t75 = cos(t77);
	t86 = t80 * t75;
	t82 = cos(qJ(1));
	t85 = t82 * t73;
	t84 = t82 * t75;
	t72 = cos(qJ(4)) * pkin(3) + pkin(2);
	t78 = qJ(2) + qJ(3);
	t74 = sin(t78);
	t76 = cos(t78);
	t83 = pkin(1) * cos(qJ(2)) + pkin(5) * t74 + t72 * t76;
	t1 = [t76 * t84 + t87, -t76 * t85 + t86, t82 * t74, t80 * t88 + t83 * t82 + 0; t76 * t86 - t85, -t76 * t87 - t84, t80 * t74, t83 * t80 - t82 * t88 + 0; t74 * t75, -t74 * t73, -t76, t74 * t72 - t76 * pkin(5) + sin(qJ(2)) * pkin(1) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end