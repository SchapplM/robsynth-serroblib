% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP2
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:51
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRP2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:51:41
% EndTime: 2019-02-22 10:51:41
% DurationCPUTime: 0.04s
% Computational Cost: add. (86->18), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->17)
t89 = qJ(4) + qJ(5);
t86 = sin(t89);
t90 = sin(qJ(3));
t95 = t90 * t86;
t87 = cos(t89);
t94 = t90 * t87;
t91 = cos(qJ(3));
t93 = t91 * t86;
t92 = t91 * t87;
t88 = qJ(1) + pkin(10);
t85 = cos(t88);
t84 = sin(t88);
t83 = t84 * t86 + t85 * t92;
t82 = t84 * t87 - t85 * t93;
t81 = -t84 * t92 + t85 * t86;
t80 = t84 * t93 + t85 * t87;
t1 = [t81, 0, -t85 * t94, t82, t82, 0; t83, 0, -t84 * t94, -t80, -t80, 0; 0, 0, t92, -t95, -t95, 0; t80, 0, t85 * t95, -t83, -t83, 0; t82, 0, t84 * t95, t81, t81, 0; 0, 0, -t93, -t94, -t94, 0; -t84 * t90, 0, t85 * t91, 0, 0, 0; t85 * t90, 0, t84 * t91, 0, 0, 0; 0, 0, t90, 0, 0, 0;];
JR_rot  = t1;
