% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:43
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPPR3_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_jacobiR_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:43:18
% EndTime: 2019-02-22 09:43:19
% DurationCPUTime: 0.08s
% Computational Cost: add. (16->10), mult. (57->32), div. (0->0), fcn. (88->8), ass. (0->19)
t84 = sin(pkin(6));
t87 = sin(qJ(3));
t96 = t84 * t87;
t88 = sin(qJ(2));
t95 = t84 * t88;
t89 = cos(qJ(3));
t94 = t84 * t89;
t90 = cos(qJ(2));
t93 = t84 * t90;
t86 = cos(pkin(6));
t92 = t86 * t88;
t91 = t86 * t90;
t85 = cos(pkin(10));
t83 = sin(pkin(10));
t82 = -t83 * t92 + t85 * t90;
t81 = -t83 * t91 - t85 * t88;
t80 = t83 * t90 + t85 * t92;
t79 = -t83 * t88 + t85 * t91;
t1 = [0, t81 * t89, -t82 * t87 + t83 * t94, 0, 0, 0; 0, t79 * t89, -t80 * t87 - t85 * t94, 0, 0, 0; 0, t89 * t93, t86 * t89 - t87 * t95, 0, 0, 0; 0, t82, 0, 0, 0, 0; 0, t80, 0, 0, 0, 0; 0, t95, 0, 0, 0, 0; 0, t81 * t87, t82 * t89 + t83 * t96, 0, 0, 0; 0, t79 * t87, t80 * t89 - t85 * t96, 0, 0, 0; 0, t87 * t93, t86 * t87 + t88 * t94, 0, 0, 0;];
JR_rot  = t1;
