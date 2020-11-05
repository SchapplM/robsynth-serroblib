% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPRP3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:01
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPPRP3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:01:08
	% EndTime: 2020-11-04 22:01:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:01:08
	% EndTime: 2020-11-04 22:01:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t48 = cos(qJ(1));
	t47 = sin(qJ(1));
	t1 = [t48, -t47, 0, 0; t47, t48, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:01:08
	% EndTime: 2020-11-04 22:01:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t52 = cos(qJ(1));
	t51 = cos(qJ(2));
	t50 = sin(qJ(1));
	t49 = sin(qJ(2));
	t1 = [t52 * t51, -t52 * t49, t50, t52 * pkin(1) + t50 * pkin(7) + 0; t50 * t51, -t50 * t49, -t52, t50 * pkin(1) - t52 * pkin(7) + 0; t49, t51, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:01:08
	% EndTime: 2020-11-04 22:01:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (13->11), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->6)
	t57 = cos(qJ(1));
	t56 = cos(qJ(2));
	t55 = sin(qJ(1));
	t54 = sin(qJ(2));
	t53 = t56 * pkin(2) + t54 * qJ(3) + pkin(1);
	t1 = [t57 * t56, t55, t57 * t54, t55 * pkin(7) + t53 * t57 + 0; t55 * t56, -t57, t55 * t54, -t57 * pkin(7) + t53 * t55 + 0; t54, 0, -t56, t54 * pkin(2) - t56 * qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:01:08
	% EndTime: 2020-11-04 22:01:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->16), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->8)
	t64 = pkin(2) + pkin(3);
	t63 = cos(qJ(1));
	t62 = cos(qJ(2));
	t61 = sin(qJ(1));
	t60 = sin(qJ(2));
	t59 = pkin(7) - qJ(4);
	t58 = t60 * qJ(3) + t64 * t62 + pkin(1);
	t1 = [t63 * t60, -t63 * t62, -t61, t58 * t63 + t59 * t61 + 0; t61 * t60, -t61 * t62, t63, t58 * t61 - t59 * t63 + 0; -t62, -t60, 0, -t62 * qJ(3) + t64 * t60 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:01:09
	% EndTime: 2020-11-04 22:01:09
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (27->17), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->15)
	t69 = sin(qJ(5));
	t71 = sin(qJ(1));
	t78 = t71 * t69;
	t72 = cos(qJ(5));
	t77 = t71 * t72;
	t74 = cos(qJ(1));
	t76 = t74 * t69;
	t75 = t74 * t72;
	t73 = cos(qJ(2));
	t70 = sin(qJ(2));
	t68 = pkin(7) - qJ(4);
	t67 = qJ(3) + pkin(4);
	t66 = pkin(2) + pkin(3) + pkin(8);
	t65 = t66 * t73 + t67 * t70 + pkin(1);
	t1 = [t70 * t75 - t78, -t70 * t76 - t77, t74 * t73, t65 * t74 + t68 * t71 + 0; t70 * t77 + t76, -t70 * t78 + t75, t71 * t73, t65 * t71 - t68 * t74 + 0; -t73 * t72, t73 * t69, t70, t66 * t70 - t67 * t73 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:01:09
	% EndTime: 2020-11-04 22:01:09
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (35->19), mult. (31->22), div. (0->0), fcn. (44->6), ass. (0->15)
	t83 = sin(qJ(5));
	t85 = sin(qJ(1));
	t92 = t85 * t83;
	t86 = cos(qJ(5));
	t91 = t85 * t86;
	t88 = cos(qJ(1));
	t90 = t88 * t83;
	t89 = t88 * t86;
	t87 = cos(qJ(2));
	t84 = sin(qJ(2));
	t82 = qJ(6) + pkin(2) + pkin(3) + pkin(8);
	t81 = t86 * pkin(5) + pkin(4) + qJ(3);
	t80 = t83 * pkin(5) - pkin(7) + qJ(4);
	t79 = t81 * t84 + t82 * t87 + pkin(1);
	t1 = [t84 * t89 - t92, -t84 * t90 - t91, t88 * t87, t79 * t88 - t80 * t85 + 0; t84 * t91 + t90, -t84 * t92 + t89, t85 * t87, t79 * t85 + t80 * t88 + 0; -t87 * t86, t87 * t83, t84, -t81 * t87 + t82 * t84 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end