% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR3
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:59
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRR3_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:59:35
% EndTime: 2019-02-22 10:59:35
% DurationCPUTime: 0.04s
% Computational Cost: add. (86->18), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->17)
t88 = qJ(4) + qJ(5);
t85 = sin(t88);
t89 = sin(qJ(3));
t94 = t89 * t85;
t86 = cos(t88);
t93 = t89 * t86;
t90 = cos(qJ(3));
t92 = t90 * t85;
t91 = t90 * t86;
t87 = qJ(1) + pkin(11);
t84 = cos(t87);
t83 = sin(t87);
t82 = t83 * t85 + t84 * t91;
t81 = t83 * t86 - t84 * t92;
t80 = -t83 * t91 + t84 * t85;
t79 = t83 * t92 + t84 * t86;
t1 = [t80, 0, -t84 * t93, t81, t81, 0; t82, 0, -t83 * t93, -t79, -t79, 0; 0, 0, t91, -t94, -t94, 0; t79, 0, t84 * t94, -t82, -t82, 0; t81, 0, t83 * t94, t80, t80, 0; 0, 0, -t92, -t93, -t93, 0; -t83 * t89, 0, t84 * t90, 0, 0, 0; t84 * t89, 0, t83 * t90, 0, 0, 0; 0, 0, t89, 0, 0, 0;];
JR_rot  = t1;
