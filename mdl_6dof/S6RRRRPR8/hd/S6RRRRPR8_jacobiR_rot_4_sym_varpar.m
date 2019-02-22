% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:21
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR8_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_jacobiR_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:21:32
% EndTime: 2019-02-22 12:21:32
% DurationCPUTime: 0.04s
% Computational Cost: add. (54->17), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->19)
t80 = qJ(3) + qJ(4);
t78 = sin(t80);
t81 = sin(qJ(2));
t91 = t81 * t78;
t79 = cos(t80);
t90 = t81 * t79;
t82 = sin(qJ(1));
t89 = t82 * t81;
t83 = cos(qJ(2));
t88 = t83 * t78;
t87 = t83 * t79;
t84 = cos(qJ(1));
t86 = t84 * t81;
t85 = t84 * t83;
t77 = t82 * t78 + t79 * t85;
t76 = -t78 * t85 + t82 * t79;
t75 = t84 * t78 - t82 * t87;
t74 = t84 * t79 + t82 * t88;
t1 = [t75, -t79 * t86, t76, t76, 0, 0; t77, -t79 * t89, -t74, -t74, 0, 0; 0, t87, -t91, -t91, 0, 0; t74, t78 * t86, -t77, -t77, 0, 0; t76, t78 * t89, t75, t75, 0, 0; 0, -t88, -t90, -t90, 0, 0; -t89, t85, 0, 0, 0, 0; t86, t82 * t83, 0, 0, 0, 0; 0, t81, 0, 0, 0, 0;];
JR_rot  = t1;
