% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRRR8 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:32
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPPRRR8_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR8_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:32:44
	% EndTime: 2020-11-04 21:32:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:32:44
	% EndTime: 2020-11-04 21:32:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t54 = cos(qJ(1));
	t53 = sin(qJ(1));
	t1 = [t54, -t53, 0, 0; t53, t54, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:32:44
	% EndTime: 2020-11-04 21:32:44
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t56 = cos(qJ(1));
	t55 = sin(qJ(1));
	t1 = [0, -t56, t55, t56 * pkin(1) + t55 * qJ(2) + 0; 0, -t55, -t56, t55 * pkin(1) - t56 * qJ(2) + 0; 1, 0, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:32:44
	% EndTime: 2020-11-04 21:32:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t61 = cos(qJ(1));
	t60 = sin(qJ(1));
	t59 = pkin(1) + qJ(3);
	t58 = cos(pkin(10));
	t57 = sin(pkin(10));
	t1 = [t60 * t57, t60 * t58, t61, t60 * qJ(2) + t59 * t61 + 0; -t61 * t57, -t61 * t58, t60, -t61 * qJ(2) + t59 * t60 + 0; t58, -t57, 0, pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:32:44
	% EndTime: 2020-11-04 21:32:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->13), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t68 = cos(qJ(1));
	t67 = sin(qJ(1));
	t66 = pkin(10) + qJ(4);
	t65 = pkin(1) + pkin(7) + qJ(3);
	t64 = cos(t66);
	t63 = sin(t66);
	t62 = sin(pkin(10)) * pkin(3) + qJ(2);
	t1 = [t67 * t63, t67 * t64, t68, t62 * t67 + t65 * t68 + 0; -t68 * t63, -t68 * t64, t67, -t62 * t68 + t65 * t67 + 0; t64, -t63, 0, cos(pkin(10)) * pkin(3) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:32:44
	% EndTime: 2020-11-04 21:32:44
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (40->20), mult. (33->22), div. (0->0), fcn. (46->8), ass. (0->14)
	t74 = sin(qJ(5));
	t75 = sin(qJ(1));
	t82 = t75 * t74;
	t76 = cos(qJ(5));
	t81 = t75 * t76;
	t77 = cos(qJ(1));
	t80 = t77 * t74;
	t79 = t77 * t76;
	t73 = pkin(10) + qJ(4);
	t70 = sin(t73);
	t71 = cos(t73);
	t78 = pkin(4) * t70 - pkin(8) * t71 + sin(pkin(10)) * pkin(3) + qJ(2);
	t72 = pkin(1) + pkin(7) + qJ(3);
	t1 = [t70 * t81 + t80, -t70 * t82 + t79, -t75 * t71, t72 * t77 + t78 * t75 + 0; -t70 * t79 + t82, t70 * t80 + t81, t77 * t71, t72 * t75 - t78 * t77 + 0; t71 * t76, -t71 * t74, t70, t71 * pkin(4) + t70 * pkin(8) + cos(pkin(10)) * pkin(3) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:32:44
	% EndTime: 2020-11-04 21:32:44
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (58->25), mult. (40->24), div. (0->0), fcn. (53->10), ass. (0->17)
	t91 = qJ(5) + qJ(6);
	t87 = sin(t91);
	t93 = sin(qJ(1));
	t101 = t93 * t87;
	t88 = cos(t91);
	t100 = t93 * t88;
	t94 = cos(qJ(1));
	t99 = t94 * t87;
	t98 = t94 * t88;
	t97 = pkin(5) * sin(qJ(5)) + pkin(1) + pkin(7) + qJ(3);
	t84 = cos(qJ(5)) * pkin(5) + pkin(4);
	t90 = pkin(10) + qJ(4);
	t85 = sin(t90);
	t86 = cos(t90);
	t95 = -pkin(9) - pkin(8);
	t96 = t84 * t85 + t86 * t95 + sin(pkin(10)) * pkin(3) + qJ(2);
	t1 = [t85 * t100 + t99, -t85 * t101 + t98, -t93 * t86, t96 * t93 + t97 * t94 + 0; -t85 * t98 + t101, t85 * t99 + t100, t94 * t86, t97 * t93 - t96 * t94 + 0; t86 * t88, -t86 * t87, t85, t86 * t84 - t85 * t95 + cos(pkin(10)) * pkin(3) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end